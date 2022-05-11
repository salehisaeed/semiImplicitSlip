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

#include "displacementDirichletBCPointPatchVectorField.H"
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

displacementDirichletBCPointPatchVectorField::
displacementDirichletBCPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    p0_(p.localPoints())
{}


displacementDirichletBCPointPatchVectorField::
displacementDirichletBCPointPatchVectorField
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


displacementDirichletBCPointPatchVectorField::
displacementDirichletBCPointPatchVectorField
(
    const displacementDirichletBCPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    p0_(ptf.p0_, mapper)
{}


displacementDirichletBCPointPatchVectorField::
displacementDirichletBCPointPatchVectorField
(
    const displacementDirichletBCPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void displacementDirichletBCPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}


void displacementDirichletBCPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const displacementDirichletBCPointPatchVectorField& aODptf =
        refCast<const displacementDirichletBCPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}


void displacementDirichletBCPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Read mesh points
    const polyMesh& mesh = this->internalField().mesh()();
    pointField points(mesh.points());

    // Read the pointDisplacement0_ field
    const pointVectorField& pointDisplacement0_ =
        this->db().objectRegistry::lookupObject
        <pointVectorField> ("pointDisplacement0");

    vectorField displacement(this->patchInternalField());
    labelList patchPoints = patch().meshPoints();

    forAll(displacement, idx)//loop over all patch points
    {
        displacement[idx] = pointDisplacement0_[patchPoints[idx]];
    }


    // ********************************************************************** //
    // Correct the displacements to make sure that the points are following the
    // exact mathematical profile of the hub
    scalar xCoor, yCoor;
    scalar radii, radii2, radii3;
    scalar zCoor, zCoorNew, zCoor1, zCoor2, zCoor3;

    //radii1 = 0.28;
    radii2 = 0.22985032258099999;
    radii3 = 0.20735321320750852;
    //radii4 = 0.20274838709707002;
    zCoor1 = 0.083354838709679996;
    zCoor2 = 0.083354838709840007;
    zCoor3 = 0.096011464515818892;
    //zCoor4 = 0.106064516129023540;

    int n = 2;
    std::vector<double> hrc(n,0), hzc(n,0), hCircRadii(n,0);

    hrc[0] = 0.22985032258099999;
    hrc[1] = 0.25264516129000003;

    hzc[0] = 0.109677419355000010;
    hzc[1] = 0.122838709677000040;

    hCircRadii[0] = Foam::sqrt(pow(hrc[0]-radii2,2) + pow(hzc[0]-zCoor2,2));
    hCircRadii[1] = Foam::sqrt(pow(hrc[1]-radii3,2) + pow(hzc[1]-zCoor3,2));

    forAll(p0_, idx)
    {
        xCoor = displacement[idx].x() + p0_[idx].x();
        yCoor = displacement[idx].y() + p0_[idx].y();
        zCoor = displacement[idx].z() + p0_[idx].z();

        radii = Foam::sqrt(pow(xCoor,2) + pow(yCoor,2));

        if (radii > radii2)
        {
            zCoorNew = zCoor1;
        }
        else if (radii <= radii2 && radii > radii3)
        {
            zCoorNew = hzc[0]-Foam::sqrt(pow(hCircRadii[0],2)-pow(radii-hrc[0],2));
        }
        else
        {
            zCoorNew = hzc[1]-Foam::sqrt(pow(hCircRadii[1],2)-pow(radii-hrc[1],2));
        }

        displacement[idx].z() += zCoorNew - zCoor;
    }
    // ********************************************************************** //


    vectorField::operator=
    (
        displacement
    );

    fixedValuePointPatchField<vector>::updateCoeffs();

}


void displacementDirichletBCPointPatchVectorField::write
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
    displacementDirichletBCPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
