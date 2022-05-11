/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "rotatingSlipFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::rotatingSlipFvPatchField<Type>::rotatingSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    rotatingBasicSymmetryFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::rotatingSlipFvPatchField<Type>::rotatingSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    rotatingBasicSymmetryFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::rotatingSlipFvPatchField<Type>::rotatingSlipFvPatchField
(
    const rotatingSlipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    rotatingBasicSymmetryFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::rotatingSlipFvPatchField<Type>::rotatingSlipFvPatchField
(
    const rotatingSlipFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    rotatingBasicSymmetryFvPatchField<Type>(ptf, iF)
{}


template<class Type>
Foam::rotatingSlipFvPatchField<Type>::rotatingSlipFvPatchField
(
    const rotatingSlipFvPatchField<Type>& ptf
)
:
    rotatingBasicSymmetryFvPatchField<Type>(ptf)
{}


// ************************************************************************* //
