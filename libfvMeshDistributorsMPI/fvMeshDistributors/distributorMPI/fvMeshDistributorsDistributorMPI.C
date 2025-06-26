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

#include "fvMeshDistributorsDistributorMPI.H"
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
    defineTypeNameAndDebug(distributorMPI, 0);
    addToRunTimeSelectionTable
    (
        fvMeshDistributor,
        distributorMPI,
        fvMesh
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvMeshDistributors::distributorMPI::readDict()
{
    const dictionary& distributorDict(dict());

    redistributionInterval_ =
        distributorDict.lookupOrDefault("redistributionInterval", 10);

    maxImbalance_ =
        distributorDict.lookupOrDefault<scalar>("maxImbalance", 0.1);

    maxMemImbalance_ =
        distributorDict.lookupOrDefault<scalar>("maxMemImbalance", 0.1);
}


void Foam::fvMeshDistributors::distributorMPI::distribute
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

Foam::fvMeshDistributors::distributorMPI::distributorMPI(fvMesh& mesh)
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
    timeIndex_(-1)
{
    readDict();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshDistributors::distributorMPI::~distributorMPI()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshDistributors::distributorMPI::update()
{
    const fvMesh& mesh = this->mesh();

    bool redistributed = false;

    if
    (
        Pstream::nProcs() > 1
     && mesh.time().timeIndex() > 1
     && timeIndex_ != mesh.time().timeIndex()
     && mesh.time().timeIndex() % redistributionInterval_ == 0
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


        if (loadImbalance_ > maxImbalance_ || memoryImbalance_ > maxMemImbalance_)
        {
            Info<< "Redistributing mesh due to imbalance!" << endl;
            Info<< "Max CPU imbalance: " << loadImbalance_ << endl;
            Info<< "Max MEM imbalance: " << memoryImbalance_ << endl;

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


void Foam::fvMeshDistributors::distributorMPI::topoChange(const polyTopoChangeMap&)
{}


void Foam::fvMeshDistributors::distributorMPI::mapMesh(const polyMeshMap&)
{}


void Foam::fvMeshDistributors::distributorMPI::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fvMeshDistributors::distributorMPI::write(const bool write) const
{
    return true;
}


// ************************************************************************* //
