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

\*---------------------------------------------------------------------------*/

#include "rotatingBasicSymmetryPointPatchField.H"
#include "transformField.H"
#include "symmTransformField.H"
#include "mathematicalConstants.H"
#include "Time.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::rotatingBasicSymmetryPointPatchField<Type>::rotatingBasicSymmetryPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    pointPatchField<Type>(p, iF),
    origin_(),
    axis_(Zero),    
    omega_(0.0)
{}


template<class Type>
Foam::rotatingBasicSymmetryPointPatchField<Type>::rotatingBasicSymmetryPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    pointPatchField<Type>(p, iF, dict),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),    
    omega_(dict.get<scalar>("omega"))
{}


template<class Type>
Foam::rotatingBasicSymmetryPointPatchField<Type>::rotatingBasicSymmetryPointPatchField
(
    const rotatingBasicSymmetryPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    pointPatchField<Type>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    omega_(ptf.omega_)
{}


template<class Type>
Foam::rotatingBasicSymmetryPointPatchField<Type>::rotatingBasicSymmetryPointPatchField
(
    const rotatingBasicSymmetryPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    pointPatchField<Type>(ptf, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),    
    omega_(ptf.omega_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::rotatingBasicSymmetryPointPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    const vectorField& nHatRotated = this->patch().pointNormals();

    const polyMesh& mesh = this->internalField().mesh()();
    scalar t = mesh.time().value();
    scalar dt = mesh.time().deltaTValue();

    scalar angle = -omega_*(t - dt);
    vector axisHat = axis_/mag(axis_);
    vectorField nHatRotatedRel(nHatRotated - origin_);

    vectorField nHatOriginal
    (
        nHatRotatedRel*cos(angle)
        + (axisHat ^ nHatRotatedRel*sin(angle))
        + (axisHat & nHatRotatedRel)*(1 - cos(angle))*axisHat
    );

    nHatOriginal = nHatOriginal / mag(nHatOriginal);

    tmp<Field<Type>> tvalues =
    (
        (
            this->patchInternalField()
          + transform(I - 2.0*sqr(nHatOriginal), this->patchInternalField())
        )/2.0
    );

    // Get internal field to insert values into
    Field<Type>& iF = const_cast<Field<Type>&>(this->primitiveField());

    this->setInInternalField(iF, tvalues());
}

template<class Type>
void Foam::rotatingBasicSymmetryPointPatchField<Type>::write
(
    Ostream& os
) const
{
    pointPatchField<Type>::write(os);
    os.writeEntry("axis", axis_);
    os.writeEntry("origin", origin_);    
    os.writeEntry("omega", omega_);
}


// ************************************************************************* //
