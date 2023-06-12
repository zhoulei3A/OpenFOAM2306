/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::abaqusMeshSet

Description
    For a dictionary specification:
    \table
        Property | Description                             | Required | Default
        type     | abaqusMesh                              | yes      |
        scale    | scale                                   | no       | 1
        file     | Path to Abaqus file                     | yes      |
    \endtable

SourceFiles
    abaqusMeshSet.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_abaqusMeshSet_H
#define Foam_abaqusMeshSet_H

#include "sampledSet.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class abaqusMeshSet Declaration
\*---------------------------------------------------------------------------*/

class abaqusMeshSet
:
    public sampledSet
{
    // Private data

        //- Scale
        const scalar scale_;

        //- Sample coordinates
        vectorList sampleCoords_;

        //- Maximum search distance-squared
        const scalar maxDistSqr_;


    // Private Member Functions

        bool readCoord(ISstream& is, vector& coord) const;

        //- Samples all points in sampleCoords.
        void calcSamples
        (
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Uses calcSamples to obtain samples. Copies them into *this.
        void genSamples();


public:

    //- Runtime type information
    TypeName("abaqusMesh");


    // Constructors

        //- Construct from dictionary
        abaqusMeshSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~abaqusMeshSet() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
