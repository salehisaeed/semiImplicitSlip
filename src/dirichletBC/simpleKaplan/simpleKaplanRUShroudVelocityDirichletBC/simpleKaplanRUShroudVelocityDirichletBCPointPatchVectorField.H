/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField

Description
    Foam::simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField

SourceFiles
    simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField_H
#define simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

        pointField p0_;


public:

    //- Runtime type information
    TypeName("simpleKaplanRUShroudVelocityDirichletBC");


    // Constructors

        //- Construct from patch and internal field
        simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField
        (
            const simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField
        (
            const simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new simpleKaplanRUShroudVelocityDirichletBCPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            );

            //- Reverse map the given pointPatchField onto this pointPatchField
            virtual void rmap
            (
                const pointPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
