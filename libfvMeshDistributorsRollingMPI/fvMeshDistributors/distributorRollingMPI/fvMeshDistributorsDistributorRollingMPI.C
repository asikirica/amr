/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "fvMeshDistributorsDistributorRollingMPI.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "polyDistributionMap.H"
#include "addToRunTimeSelectionTable.H"

#include "UPstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshDistributors
{
    defineTypeNameAndDebug(distributorRollingMPI, 0);
    addToRunTimeSelectionTable
    (
        fvMeshDistributor,
        distributorRollingMPI,
        fvMesh
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvMeshDistributors::distributorRollingMPI::readDict()
{
    const dictionary& distributorDict(dict());

    redistributionInterval_ =
        distributorDict.lookupOrDefault("redistributionInterval", 10);

    maxImbalance_ =
        distributorDict.lookupOrDefault<scalar>("maxImbalance", 0.1);

    maxMemImbalance_ =
        distributorDict.lookupOrDefault<scalar>("maxMemImbalance", 0.1);

    maxImbalancedStates_ =
        min
        (
            redistributionInterval_,
            max
            (
                1,
                distributorDict.lookupOrDefault<label>
                (
                    "maxImbalancedStates",
                    redistributionInterval_/2
                )
            )
        );
}


void Foam::fvMeshDistributors::distributorRollingMPI::distribute
(
    const labelList& distribution
)
{
    fvMesh& mesh = this->mesh();

    // Mesh distribution engine
    fvMeshDistribute distributor(mesh);

    // Do actual sending/receiving of mesh
    autoPtr<polyDistributionMap> map
    (
        distributor.distribute(distribution)
    );

    // Distribute the mesh data
    mesh.distribute(map);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::distributorRollingMPI::distributorRollingMPI(fvMesh& mesh)
:
    fvMeshDistributor(mesh),
    distributor_
    (
        decompositionMethod::NewDistributor
        (
            decompositionMethod::decomposeParDict(mesh.time())
        )
    ),
    redistributionInterval_(1),
    maxImbalance_(0.1),
    maxMemImbalance_(0.1),
    maxImbalancedStates_(1),
    timeIndex_(-1),
    historyIndex_(0)
{
    readDict();

    loadImbalanceHistory_.setSize(redistributionInterval_, 0.0);
    memImbalanceHistory_.setSize(redistributionInterval_, 0.0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::distributorRollingMPI::~distributorRollingMPI()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshDistributors::distributorRollingMPI::update()
{
    const fvMesh& mesh = this->mesh();

    bool redistributed = false;

    if
    (
        Pstream::nProcs() > 1
     && mesh.time().timeIndex() > 1
     && timeIndex_ != mesh.time().timeIndex()
    )
    {
        timeIndex_ = mesh.time().timeIndex();

        // Get combined stats in a single call
        scalarList combinedStats = Foam::getMPIStats();


        // Get CPU load statistics
        scalarList MPIload(Pstream::nProcs());
        for(label i = 0; i < Pstream::nProcs(); i++)
        {
            MPIload[i] = combinedStats[2*i];
        }
        MPIload = returnReduce(MPIload, maxOp<scalarList>());

        // Average load imbalance
        const scalar avgMPIload = max(SMALL, average(MPIload));
        scalarList loadImbalance = mag((MPIload - avgMPIload)/avgMPIload);


        // Get memory load statistics
        scalarList MPImemory(Pstream::nProcs());
        for(label i = 0; i < Pstream::nProcs(); i++)
        {
            MPImemory[i] = combinedStats[2*i + 1];
        }
        MPImemory = returnReduce(MPImemory, maxOp<scalarList>());

        // Average memory imbalance
        const scalar avgMPImemory = max(SMALL, average(MPImemory));
        scalarList memoryImbalance = mag((MPImemory - avgMPImemory)/avgMPImemory);


        // Compute imbalance
        const scalar loadImbalance_ = max(loadImbalance);
        const scalar memoryImbalance_ = max(memoryImbalance);


        // Update history
        loadImbalanceHistory_[historyIndex_] = loadImbalance_;
        memImbalanceHistory_[historyIndex_] = memoryImbalance_;
        historyIndex_ = (historyIndex_ + 1) % redistributionInterval_;


        // Imbalanced states
        label nImbalancedStates = 0;
        for(label i=0; i < redistributionInterval_; i++)
        {
            if (loadImbalanceHistory_[i] > maxImbalance_ || memImbalanceHistory_[i] > maxMemImbalance_)
            {
                nImbalancedStates++;
            }
        }


        if (historyIndex_ == 0 && nImbalancedStates > maxImbalancedStates_)
        {
            Info<< "Redistributing mesh due to imbalance!" << endl;
            Info<< "Imbalanced states: " << nImbalancedStates << " of " << redistributionInterval_ << endl;
            Info<< "Max CPU imbalance: " << max(loadImbalanceHistory_) << endl;
            Info<< "Max MEM imbalance: " << max(memImbalanceHistory_) << endl;

            // Create new decomposition distribution
            const labelList distribution
            (
                distributor_->decompose(mesh, scalarField())
            );

            distribute(distribution);

            redistributed = true;
        }
    }

    return redistributed;
}


void Foam::fvMeshDistributors::distributorRollingMPI::topoChange(const polyTopoChangeMap&)
{}


void Foam::fvMeshDistributors::distributorRollingMPI::mapMesh(const polyMeshMap&)
{}


void Foam::fvMeshDistributors::distributorRollingMPI::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fvMeshDistributors::distributorRollingMPI::write(const bool write) const
{
    return true;
}


// ************************************************************************* //
