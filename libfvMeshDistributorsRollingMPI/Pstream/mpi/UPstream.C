/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "UPstream.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"
#include "PstreamGlobals.H"
#include "SubList.H"
#include "allReduce.H"

#include <mpi.h>
#include "memInfo.H" 

#include <cstring>
#include <cstdlib>
#include <csignal>

#include <vector>
#include <chrono>
#include <fstream>

#if defined(WM_SP)
    #define MPI_SCALAR MPI_FLOAT
#elif defined(WM_DP)
    #define MPI_SCALAR MPI_DOUBLE
#elif defined(WM_LP)
    #define MPI_SCALAR MPI_LONG_DOUBLE
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Store MPI wait statistics
    static std::vector<std::chrono::steady_clock::time_point> prevRankTimes(1);
    static std::vector<double> rankTimes(1, 0.0);
    static std::vector<double> prevRankWaitTimes(1, 0.0);
    static std::vector<double> rankWaitTimes(1, 0.0);


    // Get MPI wait statistics
    scalarList getMPIStats()
    {
        const label rank = Pstream::myProcNo();
        const label nProcs = Pstream::nProcs();
        double localStats[2];


        // Get current time
        auto current_time = std::chrono::steady_clock::now();


        // Calculate delta time
        std::chrono::duration<double> delta_time = current_time - prevRankTimes[rank];
        rankTimes[rank] = delta_time.count();
        prevRankTimes[rank] = current_time;


        // Ger current wait time
        double wait_time = (rank < static_cast<label>(rankWaitTimes.size())) 
            ? rankWaitTimes[rank] 
            : 0.0;


        // Calculate delta wait time
        double delta_wait_time = wait_time - prevRankWaitTimes[rank];
        prevRankWaitTimes[rank] = wait_time;


        // Calculate CPU load
        double cpu_load = (rankTimes[rank] - delta_wait_time) / rankTimes[rank];
        cpu_load = std::min(1.0, std::max(0.0, cpu_load));


        // Calculate MEM load
        double mem_load = static_cast<double>(memInfo().size()) / 1024.0;
  

        localStats[0] = cpu_load;
        localStats[1] = mem_load;

        scalarList allStats(2*nProcs, 0.0);

        MPI_Allgather
        (
            localStats,
            2,
            MPI_DOUBLE,
            allStats.data(),
            2,
            MPI_DOUBLE,
            MPI_COMM_WORLD
        );

        return allStats;
    }
}

// NOTE:
// valid parallel options vary between implementations, but flag common ones.
// if they are not removed by MPI_Init(), the subsequent argument processing
// will notice that they are wrong
void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
    validParOptions.insert("machinefile", "machine file");
}


bool Foam::UPstream::init(int& argc, char**& argv, const bool needsThread)
{
    // MPI_Init(&argc, &argv);
    int provided_thread_support;
    MPI_Init_thread
    (
        &argc,
        &argv,
        (
            needsThread
          ? MPI_THREAD_MULTIPLE
          : MPI_THREAD_SINGLE
        ),
        &provided_thread_support
    );

    // int numprocs;
    // MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    // int myRank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    int myGlobalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myGlobalRank);
    MPI_Comm_split
    (
        MPI_COMM_WORLD,
        1,
        myGlobalRank,
        &PstreamGlobals::MPI_COMM_FOAM
    );

    int numprocs;
    MPI_Comm_size(PstreamGlobals::MPI_COMM_FOAM, &numprocs);
    int myRank;
    MPI_Comm_rank(PstreamGlobals::MPI_COMM_FOAM, &myRank);

    if (debug)
    {
        Pout<< "UPstream::init : initialised with numProcs:" << numprocs
            << " myRank:" << myRank << endl;
    }

    if (numprocs <= 1)
    {
        FatalErrorInFunction
            << "bool IPstream::init(int& argc, char**& argv) : "
               "attempt to run parallel on 1 processor"
            << Foam::abort(FatalError);
    }


    // Initialise parallel structure
    setParRun(numprocs, provided_thread_support == MPI_THREAD_MULTIPLE);

    #ifndef SGIMPI
    string bufferSizeName = getEnv("MPI_BUFFER_SIZE");

    if (bufferSizeName.size())
    {
        int bufferSize = atoi(bufferSizeName.c_str());

        if (bufferSize)
        {
            MPI_Buffer_attach(new char[bufferSize], bufferSize);
        }
    }
    else
    {
        FatalErrorInFunction
            << "UPstream::init(int& argc, char**& argv) : "
            << "environment variable MPI_BUFFER_SIZE not defined"
            << Foam::abort(FatalError);
    }
    #endif

    // int processorNameLen;
    // char processorName[MPI_MAX_PROCESSOR_NAME];
    //
    // MPI_Get_processor_name(processorName, &processorNameLen);
    // processorName[processorNameLen] = '\0';
    // Pout<< "Processor name:" << processorName << endl;

    // Initialise rank times
    prevRankTimes.resize(numprocs, std::chrono::steady_clock::now());
    rankTimes.resize(numprocs, 0.0);
    prevRankWaitTimes.resize(numprocs, 0.0);
    rankWaitTimes.resize(numprocs, 0.0);

    return true;
}


void Foam::UPstream::exit(int errnum)
{
    if (debug)
    {
        Pout<< "UPstream::exit." << endl;
    }

    #ifndef SGIMPI
    int size;
    char* buff;
    MPI_Buffer_detach(&buff, &size);
    delete[] buff;
    #endif

    if (PstreamGlobals::outstandingRequests_.size())
    {
        label n = PstreamGlobals::outstandingRequests_.size();
        PstreamGlobals::outstandingRequests_.clear();

        WarningInFunction
            << "There are still " << n << " outstanding MPI_Requests." << endl
            << "This means that your code exited before doing a"
            << " UPstream::waitRequests()." << endl
            << "This should not happen for a normal code exit."
            << endl;
    }

    // Clean mpi communicators
    forAll(myProcNo_, communicator)
    {
        if (myProcNo_[communicator] != -1)
        {
            freePstreamCommunicator(communicator);
        }
    }

    if (errnum == 0)
    {
        MPI_Finalize();
        ::exit(errnum);
    }
    else
    {
        MPI_Abort(PstreamGlobals::MPI_COMM_FOAM, errnum);
    }
}


void Foam::UPstream::abort()
{
    MPI_Abort(PstreamGlobals::MPI_COMM_FOAM, 1);
}


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    // Get MPI rank
    const label rank = Pstream::myProcNo();

    // Start wait timer
    auto wait_start = std::chrono::steady_clock::now();

    allReduce(Value, 1, MPI_SCALAR, MPI_SUM, bop, tag, communicator);

    // End wait timer
    auto wait_end = std::chrono::steady_clock::now();

    // Calculate wait time
    std::chrono::duration<double> wait_time = wait_end - wait_start;

    // Get rank wait stats
    rankWaitTimes[rank] += wait_time.count();
}


void Foam::reduce
(
    scalar& Value,
    const minOp<scalar>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    // Get MPI rank
    const label rank = Pstream::myProcNo();

    // Start wait timer
    auto wait_start = std::chrono::steady_clock::now();

    allReduce(Value, 1, MPI_SCALAR, MPI_MIN, bop, tag, communicator);

    // End wait timer
    auto wait_end = std::chrono::steady_clock::now();

    // Calculate wait time
    std::chrono::duration<double> wait_time = wait_end - wait_start;

    // Get rank wait stats
    rankWaitTimes[rank] += wait_time.count();
}


void Foam::reduce
(
    vector2D& Value,
    const sumOp<vector2D>& bop,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }

    // Get MPI rank
    const label rank = Pstream::myProcNo();

    // Start wait timer
    auto wait_start = std::chrono::steady_clock::now();

    allReduce(Value, 2, MPI_SCALAR, MPI_SUM, bop, tag, communicator);

    // End wait timer
    auto wait_end = std::chrono::steady_clock::now();

    // Calculate wait time
    std::chrono::duration<double> wait_time = wait_end - wait_start;

    // Get rank wait stats
    rankWaitTimes[rank] += wait_time.count();
}


void Foam::sumReduce
(
    scalar& Value,
    label& Count,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    vector2D twoScalars(Value, scalar(Count));

    // Get MPI rank
    const label rank = Pstream::myProcNo();

    // Start wait timer
    auto wait_start = std::chrono::steady_clock::now();

    reduce(twoScalars, sumOp<vector2D>(), tag, communicator);

    // End wait timer
    auto wait_end = std::chrono::steady_clock::now();

    // Calculate wait time
    std::chrono::duration<double> wait_time = wait_end - wait_start;

    // Get rank wait stats
    rankWaitTimes[rank] += wait_time.count();

    Value = twoScalars.x();
    Count = twoScalars.y();
}


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator,
    label& requestID
)
{
    // Get MPI rank
    const label rank = Pstream::myProcNo();

    // Start wait timer
    auto wait_start = std::chrono::steady_clock::now();

#ifdef MPIX_COMM_TYPE_SHARED
    // Assume mpich2 with non-blocking collectives extensions. Once mpi3
    // is available this will change.
    MPI_Request request;
    scalar v = Value;
    MPIX_Ireduce
    (
        &v,
        &Value,
        1,
        MPI_SCALAR,
        MPI_SUM,
        0,              // root
        PstreamGlobals::MPICommunicators_[communicator],
        &request
    );

    requestID = PstreamGlobals::outstandingRequests_.size();
    PstreamGlobals::outstandingRequests_.append(request);

    if (UPstream::debug)
    {
        Pout<< "UPstream::allocateRequest for non-blocking reduce"
            << " : request:" << requestID
            << endl;
    }
#else
    // Non-blocking not yet implemented in mpi
    reduce(Value, bop, tag, communicator);
    requestID = -1;
#endif

    // End wait timer
    auto wait_end = std::chrono::steady_clock::now();

    // Calculate wait time
    std::chrono::duration<double> wait_time = wait_end - wait_start;

    // Get rank wait stats
    rankWaitTimes[rank] += wait_time.count();
}


void Foam::UPstream::allToAll
(
    const labelUList& sendData,
    labelUList& recvData,
    const label communicator
)
{
    label np = nProcs(communicator);

    if (sendData.size() != np || recvData.size() != np)
    {
        FatalErrorInFunction
            << "Size of sendData " << sendData.size()
            << " or size of recvData " << recvData.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        recvData.deepCopy(sendData);
    }
    else
    {
        // Get MPI rank
        const label rank = Pstream::myProcNo();

        // Start wait timer
        auto wait_start = std::chrono::steady_clock::now();

        if
        (
            MPI_Alltoall
            (
                const_cast<label*>(sendData.begin()),
                sizeof(label),
                MPI_BYTE,
                recvData.begin(),
                sizeof(label),
                MPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoall failed for " << sendData
                << " on communicator " << communicator
                << Foam::abort(FatalError);
        }

        // End wait timer
        auto wait_end = std::chrono::steady_clock::now();

        // Calculate wait time
        std::chrono::duration<double> wait_time = wait_end - wait_start;

        // Get rank wait stats
        rankWaitTimes[rank] += wait_time.count();
    }
}


void Foam::UPstream::allToAll
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,

    const label communicator
)
{
    label np = nProcs(communicator);

    if
    (
        sendSizes.size() != np
     || sendOffsets.size() != np
     || recvSizes.size() != np
     || recvOffsets.size() != np
    )
    {
        FatalErrorInFunction
            << "Size of sendSize " << sendSizes.size()
            << ", sendOffsets " << sendOffsets.size()
            << ", recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        if (recvSizes[0] != sendSizes[0])
        {
            FatalErrorInFunction
                << "Bytes to send " << sendSizes[0]
                << " does not equal bytes to receive " << recvSizes[0]
                << Foam::abort(FatalError);
        }
        memmove(recvData, &sendData[sendOffsets[0]], recvSizes[0]);
    }
    else
    {
        // Get MPI rank
        const label rank = Pstream::myProcNo();

        // Start wait timer
        auto wait_start = std::chrono::steady_clock::now();

        if
        (
            MPI_Alltoallv
            (
                const_cast<char*>(sendData),
                const_cast<int*>(sendSizes.begin()),
                const_cast<int*>(sendOffsets.begin()),
                MPI_BYTE,
                recvData,
                const_cast<int*>(recvSizes.begin()),
                const_cast<int*>(recvOffsets.begin()),
                MPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoallv failed for sendSizes " << sendSizes
                << " recvSizes " << recvSizes
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        // End wait timer
        auto wait_end = std::chrono::steady_clock::now();

        // Calculate wait time
        std::chrono::duration<double> wait_time = wait_end - wait_start;

        // Get rank wait stats
        rankWaitTimes[rank] += wait_time.count();
    }
}


void Foam::UPstream::gather
(
    const char* sendData,
    int sendSize,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,
    const label communicator
)
{
    label np = nProcs(communicator);

    if
    (
        UPstream::master(communicator)
     && (recvSizes.size() != np || recvOffsets.size() < np)
    )
    {
        // Note: allow recvOffsets to be e.g. 1 larger than np so we
        // can easily loop over the result

        FatalErrorInFunction
            << "Size of recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        memmove(recvData, sendData, sendSize);
    }
    else
    {
        // Get MPI rank
        const label rank = Pstream::myProcNo();

        // Start wait timer
        auto wait_start = std::chrono::steady_clock::now();

        if
        (
            MPI_Gatherv
            (
                const_cast<char*>(sendData),
                sendSize,
                MPI_BYTE,
                recvData,
                const_cast<int*>(recvSizes.begin()),
                const_cast<int*>(recvOffsets.begin()),
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Gatherv failed for sendSize " << sendSize
                << " recvSizes " << recvSizes
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        // End wait timer
        auto wait_end = std::chrono::steady_clock::now();

        // Calculate wait time
        std::chrono::duration<double> wait_time = wait_end - wait_start;

        // Get rank wait stats
        rankWaitTimes[rank] += wait_time.count();
    }
}


void Foam::UPstream::scatter
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    int recvSize,
    const label communicator
)
{
    label np = nProcs(communicator);

    if
    (
        UPstream::master(communicator)
     && (sendSizes.size() != np || sendOffsets.size() != np)
    )
    {
        FatalErrorInFunction
            << "Size of sendSizes " << sendSizes.size()
            << " or sendOffsets " << sendOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        memmove(recvData, sendData, recvSize);
    }
    else
    {
        // Get MPI rank
        const label rank = Pstream::myProcNo();

        // Start wait timer
        auto wait_start = std::chrono::steady_clock::now();

        if
        (
            MPI_Scatterv
            (
                const_cast<char*>(sendData),
                const_cast<int*>(sendSizes.begin()),
                const_cast<int*>(sendOffsets.begin()),
                MPI_BYTE,
                recvData,
                recvSize,
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Scatterv failed for sendSizes " << sendSizes
                << " sendOffsets " << sendOffsets
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

        // End wait timer
        auto wait_end = std::chrono::steady_clock::now();

        // Calculate wait time
        std::chrono::duration<double> wait_time = wait_end - wait_start;

        // Get rank wait stats
        rankWaitTimes[rank] += wait_time.count();
    }
}


void Foam::UPstream::allocatePstreamCommunicator
(
    const label parentIndex,
    const label index
)
{
    if (index == PstreamGlobals::MPIGroups_.size())
    {
        // Extend storage with dummy values
        MPI_Group newGroup = MPI_GROUP_NULL;
        PstreamGlobals::MPIGroups_.append(newGroup);
        MPI_Comm newComm = MPI_COMM_NULL;
        PstreamGlobals::MPICommunicators_.append(newComm);
    }
    else if (index > PstreamGlobals::MPIGroups_.size())
    {
        FatalErrorInFunction
            << "PstreamGlobals out of sync with UPstream data. Problem."
            << Foam::exit(FatalError);
    }


    if (parentIndex == -1)
    {
        // Allocate world communicator

        if (index != UPstream::worldComm)
        {
            FatalErrorInFunction
                << "world communicator should always be index "
                << UPstream::worldComm << Foam::exit(FatalError);
        }

        PstreamGlobals::MPICommunicators_[index] =
            PstreamGlobals::MPI_COMM_FOAM;
        MPI_Comm_group
        (
            PstreamGlobals::MPI_COMM_FOAM,
            &PstreamGlobals::MPIGroups_[index]
        );
        MPI_Comm_rank
        (
            PstreamGlobals::MPICommunicators_[index],
           &myProcNo_[index]
        );

        // Set the number of processes to the actual number
        int numProcs;
        MPI_Comm_size(PstreamGlobals::MPICommunicators_[index], &numProcs);

        // procIDs_[index] = identity(numProcs);
        procIDs_[index].setSize(numProcs);
        forAll(procIDs_[index], i)
        {
            procIDs_[index][i] = i;
        }
    }
    else
    {
        // Create new group
        MPI_Group_incl
        (
            PstreamGlobals::MPIGroups_[parentIndex],
            procIDs_[index].size(),
            procIDs_[index].begin(),
           &PstreamGlobals::MPIGroups_[index]
        );

        // Create new communicator
        MPI_Comm_create
        (
            PstreamGlobals::MPICommunicators_[parentIndex],
            PstreamGlobals::MPIGroups_[index],
           &PstreamGlobals::MPICommunicators_[index]
        );

        if (PstreamGlobals::MPICommunicators_[index] == MPI_COMM_NULL)
        {
            myProcNo_[index] = -1;
        }
        else
        {
            if
            (
                MPI_Comm_rank
                (
                    PstreamGlobals::MPICommunicators_[index],
                   &myProcNo_[index]
                )
            )
            {
                FatalErrorInFunction
                    << "Problem :"
                    << " when allocating communicator at " << index
                    << " from ranks " << procIDs_[index]
                    << " of parent " << parentIndex
                    << " cannot find my own rank"
                    << Foam::exit(FatalError);
            }
        }
    }
}


void Foam::UPstream::freePstreamCommunicator(const label communicator)
{
    if (communicator != UPstream::worldComm)
    {
        if (PstreamGlobals::MPICommunicators_[communicator] != MPI_COMM_NULL)
        {
            // Free communicator. Sets communicator to MPI_COMM_NULL
            MPI_Comm_free(&PstreamGlobals::MPICommunicators_[communicator]);
        }
        if (PstreamGlobals::MPIGroups_[communicator] != MPI_GROUP_NULL)
        {
            // Free greoup. Sets group to MPI_GROUP_NULL
            MPI_Group_free(&PstreamGlobals::MPIGroups_[communicator]);
        }
    }
}


Foam::label Foam::UPstream::nRequests()
{
    return PstreamGlobals::outstandingRequests_.size();
}


void Foam::UPstream::resetRequests(const label i)
{
    if (i < PstreamGlobals::outstandingRequests_.size())
    {
        PstreamGlobals::outstandingRequests_.setSize(i);
    }
}


void Foam::UPstream::waitRequests(const label start)
{
    if (debug)
    {
        Pout<< "UPstream::waitRequests : starting wait for "
            << PstreamGlobals::outstandingRequests_.size()-start
            << " outstanding requests starting at " << start << endl;
    }

    // Get MPI rank
    const label rank = Pstream::myProcNo();
    
    // Start wait timer
    auto wait_start = std::chrono::steady_clock::now();

    if (PstreamGlobals::outstandingRequests_.size())
    {
        SubList<MPI_Request> waitRequests
        (
            PstreamGlobals::outstandingRequests_,
            PstreamGlobals::outstandingRequests_.size() - start,
            start
        );

        if
        (
            MPI_Waitall
            (
                waitRequests.size(),
                waitRequests.begin(),
                MPI_STATUSES_IGNORE
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Waitall returned with error" << Foam::endl;
        }

        resetRequests(start);
    }

    // End wait timer
    auto wait_end = std::chrono::steady_clock::now();

    // Calculate wait time
    std::chrono::duration<double> wait_time = wait_end - wait_start;
    
    // Get rank wait stats
    rankWaitTimes[rank] += wait_time.count();

    if (debug)
    {
        Pout<< "UPstream::waitRequests : finished wait." << endl;
    }
}


void Foam::UPstream::waitRequest(const label i)
{
    if (debug)
    {
        Pout<< "UPstream::waitRequest : starting wait for request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    // Get MPI rank
    const label rank = Pstream::myProcNo();
    
    // Start wait timer
    auto wait_start = std::chrono::steady_clock::now();

    if
    (
        MPI_Wait
        (
           &PstreamGlobals::outstandingRequests_[i],
            MPI_STATUS_IGNORE
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Wait returned with error" << Foam::endl;
    }

    // End wait timer
    auto wait_end = std::chrono::steady_clock::now();

    // Calculate wait time
    std::chrono::duration<double> wait_time = wait_end - wait_start;
    
    // Get rank wait stats
    rankWaitTimes[rank] += wait_time.count();

    if (debug)
    {
        Pout<< "UPstream::waitRequest : finished wait for request:" << i
            << endl;
    }
}


bool Foam::UPstream::finishedRequest(const label i)
{
    if (debug)
    {
        Pout<< "UPstream::finishedRequest : checking request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    int flag;
    MPI_Test
    (
       &PstreamGlobals::outstandingRequests_[i],
       &flag,
        MPI_STATUS_IGNORE
    );

    if (debug)
    {
        Pout<< "UPstream::finishedRequest : finished request:" << i
            << endl;
    }

    return flag != 0;
}


int Foam::UPstream::allocateTag(const char* s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        // if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = 'X';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


int Foam::UPstream::allocateTag(const word& s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        // if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = 'X';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


void Foam::UPstream::freeTag(const char* s, const int tag)
{
    if (debug)
    {
        // if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = ' ';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


void Foam::UPstream::freeTag(const word& s, const int tag)
{
    if (debug)
    {
        // if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = ' ';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


// ************************************************************************* //
