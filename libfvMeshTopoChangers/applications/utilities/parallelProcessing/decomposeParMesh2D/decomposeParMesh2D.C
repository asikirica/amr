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

Application
    decomposeParMesh2D

Description
    Decomposes a mesh.

    Decomposes refinementHistory of a case for parallel
    execution of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "fvCFD.H"
#include "regionProperties.H"

#include "refinementHistory4.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeRefinementHistory
(
    const Time& database,
    const List<refinementHistory4::splitCell4>& splitCells,
    const labelList& visibleCells
)
{
    // Write refinementHistory used in adaptive mesh refinement
    refinementHistory4
    (
        IOobject
        (
            "refinementHistory4",
            database.findInstance(polyMesh::meshSubDir, "cellProcAddressing"),
            polyMesh::meshSubDir,
            database,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        splitCells,
        visibleCells,
        true
    ).write();

    Info<< "Writing refinement history to "
        << database.caseName()/database.timeName() << endl;
}


void markParentCells
(
    const label startIndex,
    const DynamicList<refinementHistory4::splitCell4>& splitCells,
    labelList& onProcCells
)
{
    label parentIndex = splitCells[startIndex].parent_;
    
    while (parentIndex != -1 && onProcCells[parentIndex] == 0)
    {
        onProcCells[parentIndex] = 1;
        parentIndex = splitCells[parentIndex].parent_;
    }
}


int main(int argc, char *argv[])
{
    argList::addNote("decompose a mesh");
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"

    const wordList regionNames(selectRegionNames(args, runTime));
    if (regionNames.size() > 1)
    {
        Info<< "Operating on regions " << regionNames[0];
        for (label regioni = 1; regioni < regionNames.size() - 1; ++ regioni)
        {
            Info<< ", " << regionNames[regioni];
        }
        Info<< " and " << regionNames.last() << nl << endl;
    }
    else if (regionNames[0] != polyMesh::defaultRegion)
    {
        Info<< "Operating on region " << regionNames[0] << nl << endl;
    }

    label nProcs = fileHandler().nProcs(args.path());
    Info<< "Found " << nProcs << " processor directories" << nl << endl;

    // Read all time databases
    PtrList<Time> databases(nProcs);
    forAll(databases, proci)
    {
        Info<< "Reading database "
            << args.caseName()/fileName(word("processor") + name(proci))
            << endl;

        databases.set
        (
            proci,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/fileName(word("processor") + name(proci))
            )
        );
    }

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];

        Info<< "\nDecomposing mesh " << regionName << nl << endl;

        // Load mesh
        fvMesh masterMesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            )
        );

        // Check for refinementHistory
        if (!isFile(runTime.findInstance(polyMesh::meshSubDir, "refinementHistory4")/polyMesh::meshSubDir/"refinementHistory4"))
        {
            FatalErrorIn("refinementHistory")
                << "Cannot find refinementHistory"
                << exit(FatalError);
        }

        // Load refinementHistory used in adaptive mesh refinement
        refinementHistory4 data
        (
            IOobject
            (
                "refinementHistory4",
                runTime.timeName(),
                polyMesh::meshSubDir,
                masterMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );


        forAll(databases, proci)
        {
            labelIOList cellProcAddressing
            (
                IOobject
                (
                    "cellProcAddressing",
                    databases[proci].findInstance(polyMesh::meshSubDir, "cellProcAddressing"),
                    polyMesh::meshSubDir,
                    databases[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

            // Storage for visibleCells
            labelList visibleCells(cellProcAddressing.size(), -1);

            // Storage for onProcCells
            labelList onProcCells(data.splitCells().size(), 0);

            // Mark visible cells
            forAll(cellProcAddressing, i)
            {
                label splitIndex = data.visibleCells()[cellProcAddressing[i]];
                visibleCells[i] = splitIndex;
                
                if (splitIndex != -1)
                {
                    onProcCells[splitIndex] = 1;
                    markParentCells
                    (
                        splitIndex,
                        data.splitCells(),
                        onProcCells
                    );
                }
            }

            // Count removed cells and adjust indices
            labelList nRemovedCells(data.splitCells().size(), 0);
            label removedCellCount = 0;

            // Calculate cumulative removed cells
            forAll(nRemovedCells, i)
            {
                if (onProcCells[i] == 0)
                {
                    removedCellCount++;
                }
                nRemovedCells[i] = removedCellCount;
            }

            // Calculate final number of split cells
            label nSplitCells = data.splitCells().size() - removedCellCount;

            // Adjust visible cell indices
            forAll(visibleCells, iCell)
            {
                if (visibleCells[iCell] != -1)
                {
                    visibleCells[iCell] -= nRemovedCells[visibleCells[iCell]];
                }
            }

            // Split cells for this processor
            refinementHistory4::splitCell4 empty;
            List<refinementHistory4::splitCell4> splitCellsProc(nSplitCells, empty);

            // Process cells that belong to this processor
            label destIndex = 0;
            forAll(data.splitCells(), srcIndex)
            {
                if (onProcCells[srcIndex] != 1) continue;

                // Copy cell data
                auto& destCell = splitCellsProc[destIndex];
                destCell = data.splitCells()[srcIndex];
                
                // Adjust parent index if exists
                if (destCell.parent_ != -1)
                {
                    destCell.parent_ -= nRemovedCells[destCell.parent_];
                }
                
                // Adjust child cell indices if they exist
                if (destCell.addedCellsPtr_.valid())
                {
                    auto& childCells = destCell.addedCellsPtr_();

                    // Update valid child indices
                    forAll(childCells, childIdx)
                    {
                        if (childCells[childIdx] >= 0)
                        {
                            childCells[childIdx] -= nRemovedCells[childCells[childIdx]];
                        }
                    }
                }
                destIndex++;
            }

            // Write refinementHistory
            writeRefinementHistory
            (
                databases[proci],
                splitCellsProc,
                visibleCells
            );
        }
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
