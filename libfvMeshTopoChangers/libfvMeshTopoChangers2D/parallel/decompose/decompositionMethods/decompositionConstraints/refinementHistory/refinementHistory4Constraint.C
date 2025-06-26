/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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

#include "refinementHistory4Constraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "refinementHistory4.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(refinementHistory4Constraint);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        refinementHistory4Constraint,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementHistory4Constraint::refinementHistory4Constraint
(
    const dictionary& constraintsDict,
    const word& modelType
)
:
    decompositionConstraint(constraintsDict, typeName)
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : setting constraints to preserve refinement history"
            << endl;
    }
}


Foam::refinementHistory4Constraint::refinementHistory4Constraint()
:
    decompositionConstraint(dictionary(), typeName)
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : setting constraints to refinement history"
            << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::refinementHistory4Constraint::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    autoPtr<const refinementHistory4> storagePtr;
    refinementHistory4 const* refPtr = nullptr;

    if (mesh.foundObject<refinementHistory4>("refinementHistory4"))
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " add : found refinementHistory4" << endl;
        }
        refPtr = &mesh.lookupObject<refinementHistory4>("refinementHistory4");
    }
    else
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " add : reading refinementHistory4 from time "
                << mesh.facesInstance() << endl;
        }
        storagePtr.reset
        (
            new refinementHistory4
            (
                IOobject
                (
                    "refinementHistory4",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh.nCells()
            )
        );
    }

    const refinementHistory4& history =
    (
        storagePtr.valid()
      ? storagePtr()
      : *refPtr
    );

    if (history.active())
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " add : adding refinementHistory4 " << endl;
        }

        // refinementHistory4 itself implements decompositionConstraint
        history.add
        (
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            explicitConnections
        );
    }
}


void Foam::refinementHistory4Constraint::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    autoPtr<const refinementHistory4> storagePtr;
    refinementHistory4 const* refPtr = nullptr;

    if (mesh.foundObject<refinementHistory4>("refinementHistory4"))
    {
        if (decompositionConstraint::debug)
        {
           Info<< type() << " apply : found refinementHistory4" << endl;
        }
        refPtr = &mesh.lookupObject<refinementHistory4>("refinementHistory4");
    }
    else
    {
        if (decompositionConstraint::debug)
        {
           Info<< type() << " apply : reading refinementHistory4 from time "
               << mesh.facesInstance() << endl;
        }
        storagePtr.reset
        (
            new refinementHistory4
            (
                IOobject
                (
                    "refinementHistory4",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh.nCells()
            )
        );
    }

    const refinementHistory4& history =
    (
        storagePtr.valid()
      ? storagePtr()
      : *refPtr
    );

    if (history.active())
    {
        if (decompositionConstraint::debug)
        {
            Info<< type() << " apply : adding refinementHistory4 " << endl;
        }

        // refinementHistory4 itself implements decompositionConstraint
        history.apply
        (
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            explicitConnections,
            decomposition
        );
    }
}


// ************************************************************************* //
