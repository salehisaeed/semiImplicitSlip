/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
semiImplicitSlipVelocityLaplacianSolidBody
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "solidBodyMotionFunction.H"
#include "fvOptions.H"
#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver::semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    velocityMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    SBMFPtr_(solidBodyMotionFunction::New(coeffDict(), mesh.time())),
    points0_(points0MotionSolver::points0IO(mesh)),
    pointMotionU0_
    (
        IOobject
        (
            "pointMotionU0",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh)
    ),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointMotionU_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointMotionU_.boundaryField())
    ),
    cellMotionU0_
    (
        IOobject
        (
            "cellMotionU0",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointMotionU0_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointMotionU0_.boundaryField())
    ),
    pointDisplacement_
    (
        IOobject
        (
            "pointDisplacement",
            time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh),
        dimensionedVector(dimLength, Zero)
    ),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    diffusivity0Ptr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity0"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver::~semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellMotionU0_,
        pointMotionU0_
    );

    interpolationPtr_->interpolate
    (
        cellMotionU_,
        pointMotionU_
    );

    pointDisplacement_.primitiveFieldRef() +=
        fvMesh_.time().deltaTValue()*pointMotionU_.primitiveField();

    tmp<pointField> tmorhphedPoints
    (
        points0_ + pointDisplacement_.primitiveField()
    );

    twoDCorrectPoints(tmorhphedPoints.ref());

    return transformPoints(SBMFPtr_().transformation(), tmorhphedPoints.ref());
}


void Foam::semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivity0Ptr_->correct();
    pointMotionU0_.boundaryFieldRef().updateCoeffs();

    fv::options& fvOptions(fv::options::New(fvMesh_));

    const label nNonOrthCorr
    (
        getOrDefault<label>("nNonOrthogonalCorrectors", 1)
    );

    for (label i=0; i<nNonOrthCorr; ++i)
    {
        fvVectorMatrix UEqn0
        (
            fvm::laplacian
            (
                dimensionedScalar("viscosity", dimViscosity, 1.0)
              * diffusivity0Ptr_->operator()(),
                cellMotionU0_,
                "laplacian(diffusivity,cellMotionU)"
            )
         ==
            fvOptions(cellMotionU0_)
        );

        fvOptions.constrain(UEqn0);
        UEqn0.solveSegregatedOrCoupled(UEqn0.solverDict());
        fvOptions.correct(cellMotionU0_);
    }

    diffusivityPtr_->correct();
    pointMotionU_.boundaryFieldRef().updateCoeffs();

    for (label i=0; i<nNonOrthCorr; ++i)
    {
        fvVectorMatrix UEqn
        (
            fvm::laplacian
            (
                dimensionedScalar("viscosity", dimViscosity, 1.0)
              * diffusivityPtr_->operator()(),
                cellMotionU_,
                "laplacian(diffusivity,cellMotionU)"
            )
         ==
            fvOptions(cellMotionU_)
        );

        fvOptions.constrain(UEqn);
        UEqn.solveSegregatedOrCoupled(UEqn.solverDict());
        fvOptions.correct(cellMotionU_);
    }

    cellMotionU0_ = cellMotionU_;
    //cellMotionU0_.primitiveFieldRef() = cellMotionU_.primitiveField();
    //pointMotionU0_.primitiveFieldRef() = pointMotionU_.primitiveField();

}


void Foam::semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    velocityMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityPtr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity")
    );

    diffusivity0Ptr_.reset(nullptr);
    diffusivity0Ptr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity0")
    );
}


// ************************************************************************* //
