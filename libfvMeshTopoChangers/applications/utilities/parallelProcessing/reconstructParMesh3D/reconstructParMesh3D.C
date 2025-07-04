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
    reconstructParMeshAMR

Description
    Reconstructs a mesh.

    Writes point/face/cell procAddressing so afterwards reconstructPar can be
    used to reconstruct fields.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "IOobjectList.H"
#include "labelIOList.H"
#include "processorPolyPatch.H"
#include "mapAddedPolyMesh.H"
#include "polyMeshAdder.H"
#include "faceCoupleInfo.H"
#include "fvMeshAdder.H"
#include "polyTopoChange.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "regionProperties.H"

#include "refinementHistory.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<faceCoupleInfo> determineCoupledFaces
(
    const label masterMeshProcStart,
    const label masterMeshProcEnd,
    const polyMesh& masterMesh,
    const label meshToAddProcStart,
    const label meshToAddProcEnd,
    const polyMesh& meshToAdd
)
{
    const polyBoundaryMesh& masterPatches = masterMesh.boundaryMesh();
    const polyBoundaryMesh& addPatches = meshToAdd.boundaryMesh();

    DynamicList<label> masterFaces
    (
        masterMesh.nFaces() - masterMesh.nInternalFaces()
    );
    DynamicList<label> addFaces
    (
        meshToAdd.nFaces() - meshToAdd.nInternalFaces()
    );

    for
    (
        label masterProci = masterMeshProcStart;
        masterProci < masterMeshProcEnd;
        masterProci++
    )
    {
        for
        (
            label addProci = meshToAddProcStart;
            addProci < meshToAddProcEnd;
            addProci++
        )
        {
            const word masterToAddName
            (
                "procBoundary" + name(masterProci) + "to" + name(addProci)
            );
            const word addToMasterName
            (
                "procBoundary" + name(addProci) + "to" + name(masterProci)
            );

            const label masterToAddID =
                masterPatches.findPatchID(masterToAddName);
            const label addToMasterID =
                addPatches.findPatchID(addToMasterName);

            if (masterToAddID != -1 && addToMasterID != -1)
            {
                const polyPatch& masterPp = masterPatches[masterToAddID];

                forAll(masterPp, i)
                {
                    masterFaces.append(masterPp.start() + i);
                }

                const polyPatch& addPp = addPatches[addToMasterID];

                forAll(addPp, i)
                {
                    addFaces.append(addPp.start() + i);
                }
            }

            if ((masterToAddID != -1) != (addToMasterID != -1))
            {
                const label foundProci =
                    masterToAddID != -1 ? masterProci : addProci;
                const word& foundName =
                    masterToAddID != -1 ? masterToAddName : addToMasterName;

                const label missingProci =
                    masterToAddID != -1 ? addProci : masterProci;
                const word& missingName =
                    masterToAddID != -1 ? addToMasterName : masterToAddName;

                FatalErrorInFunction
                    << "Patch " << foundName << " found on processor "
                    << foundProci << " but corresponding patch "
                    << missingName << " missing on processor "
                    << missingProci << exit(FatalError);
            }
        }
    }

    masterFaces.shrink();
    addFaces.shrink();

    return autoPtr<faceCoupleInfo>
    (
        new faceCoupleInfo
        (
            masterMesh,
            masterFaces,
            meshToAdd,
            addFaces
        )
    );
}


void writeCellDistribution
(
    Time& runTime,
    const fvMesh& masterMesh,
    const labelListList& cellProcAddressing

)
{
    // Write the decomposition as labelList for use with 'manual'
    // decomposition method.
    labelIOList cellDecomposition
    (
        IOobject
        (
            "cellDecomposition",
            masterMesh.facesInstance(),
            masterMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        masterMesh.nCells()
    );

    forAll(cellProcAddressing, proci)
    {
        const labelList& pCells = cellProcAddressing[proci];
        UIndirectList<label>(cellDecomposition, pCells) = proci;
    }

    cellDecomposition.write();

    Info<< nl << "Wrote decomposition to "
        << cellDecomposition.relativeObjectPath()
        << " for use in manual decomposition." << endl;


    // Write as volScalarField for postprocessing. Change time to 0
    // if was 'constant'
    {
        const scalar oldTime = runTime.value();
        const label oldIndex = runTime.timeIndex();
        if (runTime.timeName() == runTime.constant() && oldIndex == 0)
        {
            runTime.setTime(0, oldIndex+1);
        }

        volScalarField cellDist
        (
            IOobject
            (
                "cellDist",
                runTime.timeName(),
                masterMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            masterMesh,
            dimensionedScalar(dimless, 0),
            extrapolatedCalculatedFvPatchScalarField::typeName
        );

        forAll(cellDecomposition, celli)
        {
            cellDist[celli] = cellDecomposition[celli];
        }
        cellDist.correctBoundaryConditions();

        cellDist.write();

        Info<< nl << "Wrote decomposition as volScalarField to "
            << cellDist.name() << " for use in postprocessing."
            << endl;

        // Restore time
        runTime.setTime(oldTime, oldIndex);
    }
}


void loadAllData
(
    const Time& database,
    const label nSplitCells,
    const labelList& cellProcAddressing,
    const labelList& pointProcAddressing,
    DynamicList<refinementHistory::splitCell8>& splitCells,
    labelList& visibleCells,
    labelList& pointLevel,
    labelList& cellLevel
)
{
    // Load refinementHistory used in adaptive mesh refinement
    refinementHistory refinementData
    (
        IOobject
        (
            "refinementHistory",
            database.findInstance(polyMesh::meshSubDir, "refinementHistory"),
            polyMesh::meshSubDir,
            database,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // Process splitCells
    forAll(refinementData.splitCells(), i)
    {
        refinementHistory::splitCell8 splitCell = refinementData.splitCells()[i];
        if (splitCell.parent_ != -1)
        {
            splitCell.parent_ += nSplitCells;
        }
        if (splitCell.addedCellsPtr_.valid())
        {
            for (int j = 0; j < 4; j++)
            {
                splitCell.addedCellsPtr_()[j] += nSplitCells;
            }
        }
        splitCells.append(splitCell);
    }

    // Process visibleCells 
    forAll(refinementData.visibleCells(), i)
    {
        label visibleCell = refinementData.visibleCells()[i];
        visibleCells.append(visibleCell == -1 ? -1 : visibleCell + nSplitCells);
    }

    // Load pointLevel used in adaptive mesh refinement
    labelIOList pointLevelData
    (
        IOobject
        (
            "pointLevel",
            database.findInstance(polyMesh::meshSubDir, "pointLevel"),
            polyMesh::meshSubDir,
            database,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(pointLevelData, i)
    {
        pointLevel[pointProcAddressing[i]] = pointLevelData[i];
    }

    // Load cellLevel used in adaptive mesh refinement
    labelIOList cellLevelData
    (
        IOobject
        (
            "cellLevel", 
            database.findInstance(polyMesh::meshSubDir, "cellLevel"),
            polyMesh::meshSubDir,
            database,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(cellLevelData, i)
    {
        cellLevel[cellProcAddressing[i]] = cellLevelData[i];
    }
}


void writeRefinementHistory
(
    Time& runTime,
    const fvMesh& masterMesh,
    const List<refinementHistory::splitCell8>& splitCells,
    const labelList& visibleCells
)
{
    // Write refinementHistory used in adaptive mesh refinement
    refinementHistory
    (
        IOobject
        (
            "refinementHistory",
            runTime.timeName(),
            polyMesh::meshSubDir,
            masterMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        splitCells,
        visibleCells,
        true
    ).write();

    Info<< "Writing refinement history to "
        << runTime.caseName()/runTime.timeName() << endl;
}


void writeCellLevel
(
    Time& runTime,
    const fvMesh& masterMesh,
    const labelList& cellLevel
)
{
    // Write cellLevel used in adaptive mesh refinement
    labelIOList
    (
        IOobject
        (
            "cellLevel",
            runTime.timeName(),
            polyMesh::meshSubDir,
            masterMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        cellLevel
    ).write();

    Info<< "Writing cell level to "
        << runTime.caseName()/runTime.timeName() << endl;
}


void writePointLevel
(
    Time& runTime,
    const fvMesh& masterMesh,
    const labelList& pointLevel
)
{
    // Write pointLevel used in adaptive mesh refinement
    labelIOList
    (
        IOobject
        (
            "pointLevel",
            runTime.timeName(),
            polyMesh::meshSubDir,
            masterMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointLevel
    ).write();

    Info<< "Writing point level to "
        << runTime.caseName()/runTime.timeName() << endl;
}


int main(int argc, char *argv[])
{
    argList::addNote("reconstruct a mesh");
    timeSelector::addOptions(true, true);
    argList::noParallel();
    argList::addBoolOption
    (
        "cellDist",
        "write cell distribution as a labelList - for use with 'manual' "
        "decomposition method or as a volScalarField for post-processing."
    );

    #include "addRegionOption.H"
    #include "addAllRegionsOption.H"
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

    // Use the times list from the master processor
    // and select a subset based on the command-line options
    instantList timeDirs = timeSelector::select
    (
        databases[0].times(),
        args
    );

    // Loop over all times
    forAll(timeDirs, timeI)
    {
        // Set time for global database
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // Set time for all databases
        forAll(databases, proci)
        {
            databases[proci].setTime(timeDirs[timeI], timeI);
        }

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            const word regionDir =
                regionName == polyMesh::defaultRegion
              ? word::null
              : regionName;

            IOobject facesIO
            (
                "faces",
                databases[0].timeName(),
                regionDir/polyMesh::meshSubDir,
                databases[0],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            );


            // Problem: faceCompactIOList recognises both 'faceList' and
            //          'faceCompactList' so we cannot check the type
            if (!facesIO.headerOk())
            {
                Info<< "No mesh." << nl << endl;
                continue;
            }


            // Addressing from processor to reconstructed case
            labelListList cellProcAddressing(nProcs);
            labelListList faceProcAddressing(nProcs);
            labelListList pointProcAddressing(nProcs);

            // Internal faces on the final reconstructed mesh
            label masterInternalFaces;

            // Owner addressing on the final reconstructed mesh
            labelList masterOwner;

            {
                // Construct empty mesh.
                PtrList<fvMesh> masterMesh(nProcs);

                // Read all the meshes
                for (label proci=0; proci<nProcs; proci++)
                {
                    masterMesh.set
                    (
                        proci,
                        new fvMesh
                        (
                            IOobject
                            (
                                regionName,
                                runTime.timeName(),
                                runTime,
                                IOobject::NO_READ
                            ),
                            pointField(),
                            faceList(),
                            cellList()
                        )
                    );

                    fvMesh meshToAdd
                    (
                        IOobject
                        (
                            regionName,
                            databases[proci].timeName(),
                            databases[proci]
                        ),
                        false
                    );

                    // Initialise its addressing
                    cellProcAddressing[proci] = identity(meshToAdd.nCells());
                    faceProcAddressing[proci] = identity(meshToAdd.nFaces());
                    pointProcAddressing[proci] = identity(meshToAdd.nPoints());

                    // Find shared points/faces
                    autoPtr<faceCoupleInfo> couples = determineCoupledFaces
                    (
                        proci,
                        proci,
                        masterMesh[proci],
                        proci,
                        proci,
                        meshToAdd
                    );

                    // Add elements to mesh
                    autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
                    (
                        masterMesh[proci],
                        meshToAdd,
                        couples
                    );

                    // Added processor
                    inplaceRenumber
                    (
                        map().addedCellMap(),
                        cellProcAddressing[proci]
                    );
                    inplaceRenumber
                    (
                        map().addedFaceMap(),
                        faceProcAddressing[proci]
                    );
                    inplaceRenumber
                    (
                        map().addedPointMap(),
                        pointProcAddressing[proci]
                    );
                }

                // Merge the meshes
                for (label step=2; step<nProcs*2; step*=2)
                {
                    for (label proci=0; proci<nProcs; proci+=step)
                    {
                        label next = proci + step/2;
                        if(next >= nProcs)
                        {
                            continue;
                        }

                        Info<< "Merging mesh " << proci << " with " << next
                            << endl;

                        // Find shared points/faces
                        autoPtr<faceCoupleInfo> couples = determineCoupledFaces
                        (
                            proci,
                            next,
                            masterMesh[proci],
                            next,
                            proci+step,
                            masterMesh[next]
                        );

                        // Add elements to mesh
                        autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
                        (
                            masterMesh[proci],
                            masterMesh[next],
                            couples
                        );

                        // Processors that were already in masterMesh
                        for (label mergedI=proci; mergedI<next; mergedI++)
                        {
                            inplaceRenumber
                            (
                                map().oldCellMap(),
                                cellProcAddressing[mergedI]
                            );
                            inplaceRenumber
                            (
                                map().oldFaceMap(),
                                faceProcAddressing[mergedI]
                            );
                            inplaceRenumber
                            (
                                map().oldPointMap(),
                                pointProcAddressing[mergedI]
                            );
                        }

                        // Added processor
                        for
                        (
                            label addedI=next;
                            addedI<min(proci+step, nProcs);
                            addedI++
                        )
                        {
                            inplaceRenumber
                            (
                                map().addedCellMap(),
                                cellProcAddressing[addedI]
                            );
                            inplaceRenumber
                            (
                                map().addedFaceMap(),
                                faceProcAddressing[addedI]
                            );
                            inplaceRenumber
                            (
                                map().addedPointMap(),
                                pointProcAddressing[addedI]
                            );
                        }

                        masterMesh.set(next, nullptr);
                    }
                }

                for (label proci=0; proci<nProcs; proci++)
                {
                    Info<< "\nReading mesh to add from "
                        << databases[proci].caseName()
                        << " for time = " << databases[proci].timeName() << endl;
                }


                // Save some properties on the reconstructed mesh
                masterInternalFaces = masterMesh[0].nInternalFaces();
                masterOwner = masterMesh[0].faceOwner();


                Info<< "\nWriting merged mesh to "
                    << runTime.caseName()/runTime.timeName() << endl;

                if (!masterMesh[0].write())
                {
                    FatalErrorInFunction
                        << "Failed writing polyMesh."
                        << exit(FatalError);
                }

                if (args.optionFound("cellDist"))
                {
                    writeCellDistribution
                    (
                        runTime,
                        masterMesh[0],
                        cellProcAddressing
                    );
                }

                // Check for refinementHistory
                if (!isFile(args.rootPath()/args.caseName()/"processor0"/runTime.timeName()/polyMesh::meshSubDir/"refinementHistory"))
                {
                    FatalErrorIn("refinementHistory")
                        << "Cannot find refinementHistory"
                        << exit(FatalError);
                }

                // Check for pointLevel
                if (!isFile(args.rootPath()/args.caseName()/"processor0"/runTime.timeName()/polyMesh::meshSubDir/"pointLevel"))
                {
                    FatalErrorIn("pointLevel")
                        << "Cannot find pointLevel"
                        << exit(FatalError);
                }

                // Check for cellLevel
                if (!isFile(args.rootPath()/args.caseName()/"processor0"/runTime.timeName()/polyMesh::meshSubDir/"cellLevel"))
                {
                    FatalErrorIn("cellLevel")
                        << "Cannot find cellLevel"
                        << exit(FatalError);
                }

                // Storage for splitCells
                DynamicList<refinementHistory::splitCell8> splitCells_;

                // Storage for visibleCells
                labelList visibleCells_;

                // Storage for pointLevel
                labelList pointLevel(masterMesh[0].nPoints());

                // Storage for cellLevel
                labelList cellLevel(masterMesh[0].nCells());

                for (label proci=0; proci<nProcs; proci++)
                {
                    loadAllData
                    (
                        databases[proci],
                        splitCells_.size(),
                        cellProcAddressing[proci],
                        pointProcAddressing[proci], 
                        splitCells_,
                        visibleCells_,
                        pointLevel,
                        cellLevel
                    );
                }

                // Write refinementHistory
                writeRefinementHistory
                (
                    runTime,
                    masterMesh[0],
                    splitCells_,
                    visibleCells_
                );

                // Write pointLevel
                writePointLevel
                (
                    runTime,
                    masterMesh[0],
                    pointLevel
                );

                // Write cellLevel
                writeCellLevel
                (
                    runTime,
                    masterMesh[0],
                    cellLevel
                );
            }


            // Write the addressing

            Info<< "\nReconstructing the addressing from the processor meshes"
                << " to the newly reconstructed mesh" << nl << endl;

            forAll(databases, proci)
            {
                Info<< "Reading processor " << proci << " mesh from "
                    << databases[proci].caseName() << endl;

                polyMesh procMesh
                (
                    IOobject
                    (
                        regionName,
                        databases[proci].timeName(),
                        databases[proci]
                    )
                );


                // From processor point to reconstructed mesh point

                Info<< "Writing pointProcAddressing to "
                    << databases[proci].caseName()
                      /procMesh.facesInstance()
                      /polyMesh::meshSubDir
                    << endl;

                labelIOList
                (
                    IOobject
                    (
                        "pointProcAddressing",
                        procMesh.facesInstance(),
                        polyMesh::meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false                       // Do not register
                    ),
                    pointProcAddressing[proci]
                ).write();


                // From processor face to reconstructed mesh face

                Info<< "Writing faceProcAddressing to "
                    << databases[proci].caseName()
                      /procMesh.facesInstance()
                      /polyMesh::meshSubDir
                    << endl;

                labelIOList faceProcAddr
                (
                    IOobject
                    (
                        "faceProcAddressing",
                        procMesh.facesInstance(),
                        polyMesh::meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false                       // Do not register
                    ),
                    faceProcAddressing[proci]
                );

                // Now add turning index to faceProcAddressing.
                // See reconstructPar for meaning of turning index.
                forAll(faceProcAddr, procFacei)
                {
                    const label masterFacei = faceProcAddr[procFacei];

                    if
                    (
                       !procMesh.isInternalFace(procFacei)
                     && masterFacei < masterInternalFaces
                    )
                    {
                        // proc face is now external but used to be internal
                        // face. Check if we have owner or neighbour.

                        label procOwn = procMesh.faceOwner()[procFacei];
                        label masterOwn = masterOwner[masterFacei];

                        if (cellProcAddressing[proci][procOwn] == masterOwn)
                        {
                            // No turning. Offset by 1.
                            faceProcAddr[procFacei]++;
                        }
                        else
                        {
                            // Turned face.
                            faceProcAddr[procFacei] =
                                -1 - faceProcAddr[procFacei];
                        }
                    }
                    else
                    {
                        // No turning. Offset by 1.
                        faceProcAddr[procFacei]++;
                    }
                }

                faceProcAddr.write();


                // From processor cell to reconstructed mesh cell

                Info<< "Writing cellProcAddressing to "
                    << databases[proci].caseName()
                      /procMesh.facesInstance()
                      /polyMesh::meshSubDir
                    << endl;

                labelIOList
                (
                    IOobject
                    (
                        "cellProcAddressing",
                        procMesh.facesInstance(),
                        polyMesh::meshSubDir,
                        procMesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false                       // Do not register
                    ),
                    cellProcAddressing[proci]
                ).write();

                Info<< endl;
            }
        }
    }


    Info<< "End.\n" << endl;

    return 0;
}


// ************************************************************************* //
