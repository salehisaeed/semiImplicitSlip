/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::rotatingBasicSymmetryFvPatchField

Group
    grpConstraintBoundaryConditions

Description
    A symmetry patch

SourceFiles
    rotatingBasicSymmetryFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingBasicSymmetryFvPatchField_H
#define rotatingBasicSymmetryFvPatchField_H

#include "transformFvPatchField.H"
#include "symmetryFvPatch.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class rotatingBasicSymmetryFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class rotatingBasicSymmetryFvPatchField
:
    public transformFvPatchField<Type>
{
    // Private data

        //- Origin of the rotation
        vector origin_;

        //- Axis of the rotation
        vector axis_;

        //- Rotational speed
        scalar omega_;


public:

    //- Runtime type information
    //TypeName("rotatingBasicSymmetry");

    // Constructors

        //- Construct from patch and internal field
        rotatingBasicSymmetryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rotatingBasicSymmetryFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rotatingBasicSymmetryFvPatchField onto a new
        //  patch
        rotatingBasicSymmetryFvPatchField
        (
            const rotatingBasicSymmetryFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        rotatingBasicSymmetryFvPatchField
        (
            const rotatingBasicSymmetryFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new rotatingBasicSymmetryFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        rotatingBasicSymmetryFvPatchField
        (
            const rotatingBasicSymmetryFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new rotatingBasicSymmetryFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        //- Return gradient at boundary
        virtual tmp<Field<Type>> snGrad() const;

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Return face-gradient transform diagonal
        virtual tmp<Field<Type>> snGradTransformDiag() const;

        //- Write
        virtual void write(Ostream&) const;          
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField> rotatingBasicSymmetryFvPatchField<scalar>::snGrad() const;

template<>
void rotatingBasicSymmetryFvPatchField<scalar>::evaluate
(
    const Pstream::commsTypes commsType
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rotatingBasicSymmetryFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
