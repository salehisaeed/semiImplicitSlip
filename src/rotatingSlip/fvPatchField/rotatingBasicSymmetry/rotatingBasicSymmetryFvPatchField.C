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

#include "rotatingBasicSymmetryFvPatchField.H"
#include "symmTransformField.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::rotatingBasicSymmetryFvPatchField<Type>::rotatingBasicSymmetryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(p, iF),
    origin_(),
    axis_(Zero),
    omega_(0.0)
{}


template<class Type>
Foam::rotatingBasicSymmetryFvPatchField<Type>::rotatingBasicSymmetryFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF, dict),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),    
    omega_(dict.get<scalar>("omega"))    
{
    this->evaluate();
}


template<class Type>
Foam::rotatingBasicSymmetryFvPatchField<Type>::rotatingBasicSymmetryFvPatchField
(
    const rotatingBasicSymmetryFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),    
    omega_(ptf.omega_)
{}


template<class Type>
Foam::rotatingBasicSymmetryFvPatchField<Type>::rotatingBasicSymmetryFvPatchField
(
    const rotatingBasicSymmetryFvPatchField<Type>& ptf
)
:
    transformFvPatchField<Type>(ptf),
    origin_(ptf.origin_),
    axis_(ptf.axis_),    
    omega_(ptf.omega_)    
{}


template<class Type>
Foam::rotatingBasicSymmetryFvPatchField<Type>::rotatingBasicSymmetryFvPatchField
(
    const rotatingBasicSymmetryFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(ptf, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),    
    omega_(ptf.omega_)    
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::rotatingBasicSymmetryFvPatchField<Type>::snGrad() const
{
    const vectorField nHatRotated(this->patch().nf());

    const polyMesh& mesh = this->internalField().mesh();
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

    const Field<Type> iF(this->patchInternalField());

    return
        (transform(I - 2.0*sqr(nHatOriginal), iF) - iF)
       *(this->patch().deltaCoeffs()/2.0);
}


template<class Type>
void Foam::rotatingBasicSymmetryFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const vectorField nHatRotated(this->patch().nf());

    const polyMesh& mesh = this->internalField().mesh();
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

    const Field<Type> iF(this->patchInternalField());

    Field<Type>::operator=
    (
        (iF + transform(I - 2.0*sqr(nHatOriginal), iF))/2.0
    );

    transformFvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::rotatingBasicSymmetryFvPatchField<Type>::snGradTransformDiag() const
{
    const vectorField nHatRotated(this->patch().nf());

    const polyMesh& mesh = this->internalField().mesh();
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
    vectorField diag(nHatOriginal.size());

    diag.replace(vector::X, mag(nHatOriginal.component(vector::X)));
    diag.replace(vector::Y, mag(nHatOriginal.component(vector::Y)));
    diag.replace(vector::Z, mag(nHatOriginal.component(vector::Z)));

    return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}

template<class Type>
void Foam::rotatingBasicSymmetryFvPatchField<Type>::write
(
    Ostream& os
) const
{
    transformFvPatchField<Type>::write(os);
    os.writeEntry("axis", axis_);
    os.writeEntry("origin", origin_);    
    os.writeEntry("omega", omega_);
}

// ************************************************************************* //
