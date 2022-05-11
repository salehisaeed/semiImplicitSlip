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

#include "bumpDirichletBCPointPatchVectorField.H"
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

bumpDirichletBCPointPatchVectorField::
bumpDirichletBCPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    p0_(p.localPoints())
{}


bumpDirichletBCPointPatchVectorField::
bumpDirichletBCPointPatchVectorField
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


bumpDirichletBCPointPatchVectorField::
bumpDirichletBCPointPatchVectorField
(
    const bumpDirichletBCPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    p0_(ptf.p0_, mapper)
{}


bumpDirichletBCPointPatchVectorField::
bumpDirichletBCPointPatchVectorField
(
    const bumpDirichletBCPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void bumpDirichletBCPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}


void bumpDirichletBCPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const bumpDirichletBCPointPatchVectorField& aODptf =
        refCast<const bumpDirichletBCPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}


void bumpDirichletBCPointPatchVectorField::updateCoeffs()
{
    // if (this->updated())
    // {
    //     return;
    // }

    // Read the pointDisplacement0_ field
    const pointVectorField& pointDisplacement0_ =
        this->db().objectRegistry::
        lookupObject<pointVectorField> ("pointDisplacement0");

    vectorField displacement(this->patchInternalField());
    labelList patchPoints = patch().meshPoints();

    forAll(displacement, idx)//loop over all patch points
    {
        displacement[idx] = pointDisplacement0_[patchPoints[idx]];
    }


    // // ********************************************************************** //
    // // Correct the displacements to make sure that the points are following the
    // // exact mathematical profile of the bump
    // scalar pi = constant::mathematical::pi;
    // scalar xMinBump = 0.1;
    // scalar xMaxBump = 1.1;
    // scalar xCoor, yCoor, yCoorNew;
    //
    // forAll(p0_, idx)
    // {
    //     xCoor = displacement[idx].x() + p0_[idx].x();
    //     yCoor = displacement[idx].y() + p0_[idx].y();
    //
    //     if (xCoor <= xMinBump || xCoor >= xMaxBump)
    //     {
    //         yCoorNew = 0;
    //     }
    //     else
    //     {
    //         yCoorNew = 0.1*pow4(sin(pi*(xCoor-xMinBump)/1.0));
    //     }
    //     displacement[idx].y() += yCoorNew - yCoor;
    // }
    // // ********************************************************************** //


    vectorField::operator=
    (
        displacement
    );

    fixedValuePointPatchField<vector>::updateCoeffs();

}


void bumpDirichletBCPointPatchVectorField::write
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
    bumpDirichletBCPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
