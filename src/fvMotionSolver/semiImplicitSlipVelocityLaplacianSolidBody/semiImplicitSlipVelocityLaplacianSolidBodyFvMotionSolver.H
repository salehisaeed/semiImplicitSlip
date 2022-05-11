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
    Foam::semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian for the motion velocity.

SourceFiles
    semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver_H
#define semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver_H

#include "velocityMotionSolver.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;
class solidBodyMotionFunction;

/*---------------------------------------------------------------------------*\
              Class semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver
:
    public velocityMotionSolver,
    public fvMotionSolver
{
    // Private data

        //- Motion function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        //- access to intial points
        pointIOField points0_;

        //- Intermediate pointMotion field
        mutable pointVectorField pointMotionU0_;

        //- Intermediate cell-centre motion field
        mutable volVectorField cellMotionU_;

        //- Cell-centre motion field
        mutable volVectorField cellMotionU0_;

        //- Total pointDisplacement field
        mutable pointVectorField pointDisplacement_;
        
        //- Interpolation used to transfer cell displacement to the points
        autoPtr<motionInterpolation> interpolationPtr_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;
        autoPtr<motionDiffusivity> diffusivity0Ptr_;


    // Private Member Functions

        //- No copy construct
        semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver
        (
            const semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("semiImplicitSlipVelocityLaplacianSolidBody");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );


    //- Destructor
    ~semiImplicitSlipVelocityLaplacianSolidBodyFvMotionSolver();


    // Member Functions

        //- Return reference to the cell motion velocity field
        volVectorField& cellMotionU()
        {
            return cellMotionU_;
        }

        //- Return const reference to the cell motion velocity field
        const volVectorField& cellMotionU() const
        {
            return cellMotionU_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
