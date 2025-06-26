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

Class
    Foam::hexRef4Data

Description
    Various for reading/decomposing/reconstructing/distributing refinement
    data.

SourceFiles
    hexRef4Data.C

\*---------------------------------------------------------------------------*/

#ifndef hexRef4Data_H
#define hexRef4Data_H

#include "labelIOList.H"
#include "uniformDimensionedFields.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChangeMap;
class polyDistributionMap;
class refinementHistory4;
class fvMesh;

/*---------------------------------------------------------------------------*\
                         Class hexRef4Data Declaration
\*---------------------------------------------------------------------------*/

class hexRef4Data
{

private:

    // Private Data

        autoPtr<labelIOList> cellLevelPtr_;

        autoPtr<labelIOList> pointLevelPtr_;

        autoPtr<uniformDimensionedScalarField> level0EdgePtr_;

        autoPtr<refinementHistory4> refHistoryPtr_;


public:

    // Constructors

        //- Construct read. Has special provision for only some processors
        //  having the files so can be used in redistribution.
        hexRef4Data(const IOobject& io);

        //- Construct as subset
        hexRef4Data
        (
            const IOobject& io,
            const hexRef4Data&,
            const labelList& cellMap,
            const labelList& pointMap
        );

        //- Construct from multiple hexRef4Data
        hexRef4Data
        (
            const IOobject& io,
            const UPtrList<const labelList>& cellMaps,
            const UPtrList<const labelList>& pointMaps,
            const UPtrList<const hexRef4Data>&
        );

        //- Disallow default bitwise copy construction
        hexRef4Data(const hexRef4Data&) = delete;


    //- Destructor
    ~hexRef4Data();


    // Member Functions

        //- Parallel synchronise. This enforces valid objects on all processors
        //  (even if they don't have a mesh). Used by redistributePar.
        void sync(const IOobject& io);

        //- In-place update for topology changes
        void topoChange(const polyTopoChangeMap&);

        //- In-place distribute
        void distribute(const polyDistributionMap&);

        //- Write
        bool write() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const hexRef4Data&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
