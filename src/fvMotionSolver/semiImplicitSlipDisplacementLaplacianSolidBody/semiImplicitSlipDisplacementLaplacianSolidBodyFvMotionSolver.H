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
    Foam::semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian for the motion displacement.

SourceFiles
    semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver_H
#define semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;
class solidBodyMotionFunction;

/*---------------------------------------------------------------------------*\
             Class semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private data

        //- Motion function
        autoPtr<solidBodyMotionFunction> SBMFPtr_;

        //- my pointDisplacement field for slipping points on the U9 Shroud
        mutable pointVectorField pointDisplacement0_;

        //- Cell-centre motion field
        mutable volVectorField cellDisplacement_;

        //- Cell-centre motion field
        mutable volVectorField cellDisplacement0_;

        //- Optionally read point-position field. Used only for position
        //  boundary conditions.
        mutable autoPtr<pointVectorField> pointLocation_;

        //- Interpolation used to transfer cell displacement to the points
        autoPtr<motionInterpolation> interpolationPtr_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;
        autoPtr<motionDiffusivity> diffusivity0Ptr_;

        //- Frozen points (that are not on patches). -1 or points that are
        //  fixed to be at points0_ location
        label frozenPointsZone_;


    // Private Member Functions

        //- No copy construct
        semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver
        (
            const semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("semiImplicitSlipDisplacementLaplacianSolidBody");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );

        //- Construct from components
        semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointVectorField& pointDisplacement0,
            const pointIOField& points0
        );


    //- Destructor
    ~semiImplicitSlipDisplacementLaplacianSolidBodyFvMotionSolver();


    // Member Functions

        //- Return reference to the cell motion displacement field
        volVectorField& cellDisplacement()
        {
            return cellDisplacement_;
        }

        //- Return const reference to the cell motion displacement field
        const volVectorField& cellDisplacement() const
        {
            return cellDisplacement_;
        }

        //- Return reference to the diffusivity field
        motionDiffusivity& diffusivity();
        motionDiffusivity& diffusivity0();

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
