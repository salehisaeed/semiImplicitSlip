/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "semiImplicitSlipDisplacementLaplacianFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(semiImplicitSlipDisplacementLaplacianFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        semiImplicitSlipDisplacementLaplacianFvMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        semiImplicitSlipDisplacementLaplacianFvMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::semiImplicitSlipDisplacementLaplacianFvMotionSolver::semiImplicitSlipDisplacementLaplacianFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    pointDisplacement0_
    (
        IOobject
        (
            "pointDisplacement0",
            time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh)
    ),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    cellDisplacement0_
    (
        IOobject
        (
            "cellDisplacement0",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement0_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement0_.boundaryField())
    ),
    pointLocation_(nullptr),
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
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    )
{
    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (debug)
    {
        Info<< "semiImplicitSlipDisplacementLaplacianFvMotionSolver:" << nl
            << "    diffusivity       : " << diffusivityPtr_().type() << nl
            << "    frozenPoints zone : " << frozenPointsZone_ << endl;
    }


    if (io.typeHeaderOk<pointVectorField>(true))
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                io,
                pointMesh::New(fvMesh_)
            )
        );

        if (debug)
        {
            Info<< "semiImplicitSlipDisplacementLaplacianFvMotionSolver :"
                << " Read pointVectorField "
                << io.name()
                << " to be used for boundary conditions on points."
                << nl
                << "Boundary conditions:"
                << pointLocation_().boundaryField().types() << endl;
        }
    }
}


Foam::semiImplicitSlipDisplacementLaplacianFvMotionSolver::
semiImplicitSlipDisplacementLaplacianFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointVectorField& pointDisplacement0,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    pointDisplacement0_
    (
        IOobject
        (
            "pointDisplacement0",
            time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh)
    ),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    cellDisplacement0_
    (
        IOobject
        (
            "cellDisplacement0",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement0_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement0_.boundaryField())
    ),
    pointLocation_(nullptr),
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
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    )
{
    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (debug)
    {
        Info<< "semiImplicitSlipDisplacementLaplacianFvMotionSolver:" << nl
            << "    diffusivity       : " << diffusivityPtr_().type() << nl
            << "    frozenPoints zone : " << frozenPointsZone_ << endl;
    }


    if (io.typeHeaderOk<pointVectorField>(true))
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                io,
                pointMesh::New(fvMesh_)
            )
        );

        if (debug)
        {
            Info<< "semiImplicitSlipDisplacementLaplacianFvMotionSolver :"
                << " Read pointVectorField "
                << io.name()
                << " to be used for boundary conditions on points."
                << nl
                << "Boundary conditions:"
                << pointLocation_().boundaryField().types() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::semiImplicitSlipDisplacementLaplacianFvMotionSolver::
~semiImplicitSlipDisplacementLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::motionDiffusivity&
Foam::semiImplicitSlipDisplacementLaplacianFvMotionSolver::diffusivity()
{
    if (!diffusivityPtr_.valid())
    {
        diffusivityPtr_ = motionDiffusivity::New
        (
            fvMesh_,
            coeffDict().lookup("diffusivity")
        );
    }

    return *diffusivityPtr_;
}


Foam::motionDiffusivity&
Foam::semiImplicitSlipDisplacementLaplacianFvMotionSolver::diffusivity0()
{
    if (!diffusivity0Ptr_.valid())
    {
        diffusivity0Ptr_ = motionDiffusivity::New
        (
            fvMesh_,
            coeffDict().lookup("diffusivity0")
        );
    }

    return *diffusivity0Ptr_;
}



Foam::tmp<Foam::pointField>
Foam::semiImplicitSlipDisplacementLaplacianFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellDisplacement0_,
        pointDisplacement0_
    );

    interpolationPtr_->interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );

    if (pointLocation_.valid())
    {
        if (debug)
        {
            Info<< "semiImplicitSlipDisplacementLaplacianFvMotionSolver : applying "
                << " boundary conditions on " << pointLocation_().name()
                << " to new point location."
                << endl;
        }

        pointLocation_().primitiveFieldRef() =
            points0()
          + pointDisplacement_.primitiveField();

        pointLocation_().correctBoundaryConditions();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                pointLocation_()[pz[i]] = points0()[pz[i]];
            }
        }

        twoDCorrectPoints(pointLocation_().primitiveFieldRef());

        return tmp<pointField>(pointLocation_().primitiveField());
    }
    else
    {
        tmp<pointField> tcurPoints
        (
            points0() + pointDisplacement_.primitiveField()
        );
        pointField& curPoints = tcurPoints.ref();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                curPoints[pz[i]] = points0()[pz[i]];
            }
        }

        twoDCorrectPoints(curPoints);

        return tcurPoints;
    }
}


void Foam::semiImplicitSlipDisplacementLaplacianFvMotionSolver::solve()
{
    movePoints(fvMesh_.points());

    diffusivity0().correct();
    pointDisplacement0_.boundaryFieldRef().updateCoeffs();
    fv::options& fvOptions(fv::options::New(fvMesh_));

    fvVectorMatrix TEqn0
    (
        fvm::laplacian
        (
            dimensionedScalar("viscosity", dimViscosity, 1.0)
           *diffusivity0().operator()(),
            cellDisplacement0_,
            "laplacian(diffusivity,cellDisplacement)"
        )
     ==
        fvOptions(cellDisplacement0_)
    );

    fvOptions.constrain(TEqn0);

    //TEqn0.relax();
    TEqn0.solveSegregatedOrCoupled(TEqn0.solverDict());

    fvOptions.correct(cellDisplacement0_);

    interpolationPtr_->interpolate
    (
        cellDisplacement0_,
        pointDisplacement0_
    );
    pointDisplacement0_.correctBoundaryConditions();

    // cellDisplacement0_.relax();


    diffusivity().correct();
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    fvVectorMatrix TEqn
    (
        fvm::laplacian
        (
            dimensionedScalar("viscosity", dimViscosity, 1.0)
           *diffusivity().operator()(),
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )
     ==
        fvOptions(cellDisplacement_)
    );

    fvOptions.constrain(TEqn);

    //TEqn.relax();
    TEqn.solveSegregatedOrCoupled(TEqn.solverDict());

    fvOptions.correct(cellDisplacement_);

    interpolationPtr_->interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );
    pointDisplacement_.correctBoundaryConditions();

    //cellDisplacement_.relax();

    cellDisplacement0_ = cellDisplacement_;
    pointDisplacement0_ = pointDisplacement_;
}

/*void Foam::semiImplicitSlipDisplacementLaplacianFvMotionSolver::solve()
{
    movePoints(fvMesh_.points());

    diffusivity0().correct();
    diffusivity().correct();
    fv::options& fvOptions(fv::options::New(fvMesh_));

    pointVectorField& pointDisp = pointDisplacement0_;
    volVectorField& cellDisp = cellDisplacement0_;
    surfaceScalarField Diffusivity = diffusivity0().operator()();

    for (int i = 0; i < 2; ++i)
    {
        if (i == 1)
        {
            Diffusivity = diffusivity().operator()();
            pointDisp = pointDisplacement_;
            cellDisp = cellDisplacement_;
            Info << "Test" << endl;
        }

        pointDisp.boundaryFieldRef().updateCoeffs();

        fvVectorMatrix TEqn
        (
            fvm::laplacian
            (
                dimensionedScalar("viscosity", dimViscosity, 1.0)
               *Diffusivity,
                cellDisp,
                "laplacian(diffusivity,cellDisplacement)"
            )
         ==
            fvOptions(cellDisp)
        );

        fvOptions.constrain(TEqn);
        TEqn.solveSegregatedOrCoupled(TEqn.solverDict());
        fvOptions.correct(cellDisp);
    }
    cellDisplacement0_ = cellDisplacement_;
    pointDisplacement0_ = pointDisplacement_;
    //cellDisplacement0_.primitiveFieldRef() = cellDisplacement_.primitiveField();
    //pointDisplacement0_.primitiveFieldRef() = pointDisplacement_.primitiveField();
}*/


void Foam::semiImplicitSlipDisplacementLaplacianFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.clear();
}


// ************************************************************************* //
