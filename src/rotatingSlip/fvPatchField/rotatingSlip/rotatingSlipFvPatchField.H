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
    Foam::rotatingSlipFvPatchField

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a rotatingSlip constraint.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            rotatingSlip;
    }
    \endverbatim

SourceFiles
    rotatingSlipFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingSlipFvPatchField_H
#define rotatingSlipFvPatchField_H

#include "rotatingBasicSymmetryFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class rotatingSlipFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class rotatingSlipFvPatchField
:
    public rotatingBasicSymmetryFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("rotatingSlip");


    // Constructors

        //- Construct from patch and internal field
        rotatingSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rotatingSlipFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rotatingSlipFvPatchField onto a new patch
        rotatingSlipFvPatchField
        (
            const rotatingSlipFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        rotatingSlipFvPatchField
        (
            const rotatingSlipFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new rotatingSlipFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        rotatingSlipFvPatchField
        (
            const rotatingSlipFvPatchField<Type>&,
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
                new rotatingSlipFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Attributes

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return false;
            }


    // Member operators

        // virtual void operator=(const UList<Type>&) {}

        // virtual void operator=(const fvPatchField<Type>&) {}
        // virtual void operator+=(const fvPatchField<Type>&) {}
        // virtual void operator-=(const fvPatchField<Type>&) {}
        // virtual void operator*=(const fvPatchField<scalar>&) {}
        // virtual void operator/=(const fvPatchField<scalar>&) {}

        // virtual void operator+=(const Field<Type>&) {}
        // virtual void operator-=(const Field<Type>&) {}

        // virtual void operator*=(const Field<scalar>&) {}
        // virtual void operator/=(const Field<scalar>&) {}

        // virtual void operator=(const Type&) {}
        // virtual void operator+=(const Type&) {}
        // virtual void operator-=(const Type&) {}
        // virtual void operator*=(const scalar) {}
        // virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rotatingSlipFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
