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

#include "displacemenRotatingBCPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

displacemenRotatingBCPointPatchVectorField::
displacemenRotatingBCPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    axis_(Zero),
    origin_(Zero),
    startTime_(0.0),
    omega_(),
    printAngle_(false),
    p0_(p.localPoints())
{}


displacemenRotatingBCPointPatchVectorField::
displacemenRotatingBCPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    axis_(dict.lookup("axis")),
    origin_(dict.lookup("origin")),
    omega_(Function1<scalar>::New("omega", dict)),
    printAngle_(dict.getOrDefault<bool>("printAngle", false))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }

    if (dict.found("startTime"))
    {
        startTime_ = dict.get<scalar>("startTime");
    }
    else
    {
        startTime_ = scalar(0);
    }
}


displacemenRotatingBCPointPatchVectorField::
displacemenRotatingBCPointPatchVectorField
(
    const displacemenRotatingBCPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    startTime_(ptf.startTime_),
    omega_(ptf.omega_.clone()),
    printAngle_(ptf.printAngle_),
    p0_(ptf.p0_, mapper)
{}


displacemenRotatingBCPointPatchVectorField::
displacemenRotatingBCPointPatchVectorField
(
    const displacemenRotatingBCPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    startTime_(ptf.startTime_),
    omega_(ptf.omega_.clone()),
    printAngle_(ptf.printAngle_),
    p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void displacemenRotatingBCPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}


void displacemenRotatingBCPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const displacemenRotatingBCPointPatchVectorField& aODptf =
        refCast<const displacemenRotatingBCPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}


void displacemenRotatingBCPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    scalar angleDeg = omega_->integrate(startTime_, t.value());

	scalar angle = angleDeg * ((Foam::constant::mathematical::pi)/180.0);
    vector axisHat = axis_/mag(axis_);

    if (printAngle_)
    {
        Info << "Rotation angle of " << this->patch().name()
             << " = " << angleDeg << " deg" << endl;
    }

    vectorField p0Rel(p0_ - origin_);

    vectorField::operator=
    (
        p0Rel*(cos(angle) - 1)
      + (axisHat ^ p0Rel*sin(angle))
      + (axisHat & p0Rel)*(1 - cos(angle))*axisHat
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void displacemenRotatingBCPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("axis", axis_);
    os.writeEntry("origin", origin_);
    os.writeEntry("startTime", startTime_);
    omega_->writeData(os);
    os.writeEntry<bool>("printAngle", printAngle_);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    displacemenRotatingBCPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
