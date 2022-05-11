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

#include "simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "mathematicalConstants.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField::
simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    p0_(p.localPoints())
{}


simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField::
simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict)
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
}


simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField::
simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField
(
    const simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    p0_(ptf.p0_, mapper)
{}


simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField::
simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField
(
    const simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}


void simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField& aODptf =
        refCast<const simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}


void simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Read mesh points
    const polyMesh& mesh = this->internalField().mesh()();
    pointField points(mesh.points());
    scalar dt = mesh.time().deltaTValue();

    // Read the pointMotionU0 field
    const pointVectorField& pointMotionU0 =
        this->db().objectRegistry::
        lookupObject<pointVectorField> ("pointMotionU0");

    // Read the pointDisplacement field
    const pointVectorField& pointDisplacement =
        this->db().objectRegistry::
        lookupObject<pointVectorField> ("pointDisplacement");

    vectorField motion(this->patchInternalField());
    vectorField displacement(this->patchInternalField());
    labelList patchPoints = patch().meshPoints();

    forAll(motion, idx)//loop over all patch points
    {
        motion[idx] = pointMotionU0[patchPoints[idx]];
        displacement[idx] = pointDisplacement[patchPoints[idx]];
   }

    pointField p_ = patch().localPoints();

    // ********************************************************************** //
    // Correct the displacements to make sure that the points are following the
    // exact mathematical profile of the shroud
    scalar radii, radiiNew, radii1, rShroud;//, x2byy2, y2byx2, signx, signy;
    scalar xCoor, yCoor, zCoor, xCoorNew, yCoorNew, zCoor1, zCoor2;

    zCoor1 = -0.0427525;
    zCoor2 = 0.0427525;

    radii1 = 0.117462;
    rShroud = 0.125;

    forAll(p_, idx)
    {
        xCoor = motion[idx].x()*dt + displacement[idx].x() + p0_[idx].x();
        yCoor = motion[idx].y()*dt + displacement[idx].y() + p0_[idx].y();
        zCoor = motion[idx].z()*dt + displacement[idx].z() + p0_[idx].z();

        radii = Foam::sqrt(pow(xCoor,2) + pow(yCoor,2));

        if ((zCoor <= zCoor1) || (zCoor >= zCoor2))
        {
            radiiNew = radii1;
        }
        else
        {
            radiiNew = Foam::sqrt(pow(rShroud,2) - pow(zCoor,2));
        }

        xCoorNew = radiiNew/radii * xCoor;
        yCoorNew = radiiNew/radii * yCoor;

        motion[idx].x() += (xCoorNew - xCoor)/dt;
        motion[idx].y() += (yCoorNew - yCoor)/dt;   
    }
    // ********************************************************************** //


    vectorField::operator=
    (
        motion
    );

    fixedValuePointPatchField<vector>::updateCoeffs();

}


void simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    simpleKaplanRUHubVelocityDirichletBCPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
