# semiImplicitSlip

## General information
**semiImplicitSlip** is a mesh motion library in OpenFOAM that implements a robust semi-implicit algorithm for slipping the mesh points on curved surfaces. The algorithm includes two steps, namely, an explicit step based on the general slip condition and an implicit step based on the Dirichlet condition. It employs the Laplacian smoothing equations to spread the mesh deformation into the domain. Additionally, a solid-body rotation may be added on top of the deformed mesh, which could be useful for modelling the runner region in transient operation of Kaplan turbines which contains simultaneous mesh deformation and solid-body rotation of the mesh.

CFD simulations of numerous engineering fluid flows, such as transient operation of hydraulic turbines, involve an immensely complicated mesh motion process consisting of simultaneous mesh deformation and mesh slipping on highly curved surfaces. The available standard mesh motion methodology in OpenFOAM lacks some features to simulate this elaborate mesh motion. The current explicit slip boundary conditions in OpenFOAM fail to make large deformations of the mesh without inverting some of the cells. semiImplicitSlip library addresses this problem by developing a new dynamic mesh algorithm. The library is verified using simple, yet relevant, test cases, including a simple Kaplan turbine.

All the source codes are in the `src` folder, while the `tutorials` folder contains verification test cases. The library is tested and verified with the recent versions of [OpenFOAM](https://www.openfoam.com/), such as v1912, v2006, v2012, v2106, v2112. Please use the [Github](https://github.com/salehisaeed/semiImplicitSlip) version to receive the most recent updates.

More information about the theory, the developed codes, and the test cases are presented in the submitted manuscript.

## Installation

The library is developed based on OpenFOAM ([openfoam.com](https://www.openfoam.com/)). Therefore, a complete installation of OpenFOAM (preferably a recent version, such as OpenFOAM-v2112, is required. 

To install (compile) the semiImplicitSlip library, one needs to download source files using git:
```bash
git clone https://github.com/salehisaeed/semiImplicitSlip
```
Then, in a terminal where OpenFOAM is sourced, run the `Allwmake` script, i.e.,
```bash
cd semiImplicitSlip
./Allwmake
```
A proper compiler also needs to be available. Therefore, if the library is compiled on a cluster, one should load the corresponding compiler module before running the `Allwmake` file. To recompile the library, it is recommended to first clean the previous compilation using `Allwclean` file, i.e.,
```bash
./Allwclean
```


## Usage

The proper usage of the library is presented using different tutorial cases (inside the `tutorials` folder). Each tutorial contains a separate `README.md` file that provides a detailed description of the case.


## Files structure

### 1. src/fvMotionSolver
 
**1.1. `semiImplicitSlipDisplacementLaplacian`:** Semi-implicit slip fvMotionSolver based on Laplacian smoothing of the point displacement field

**1.2. `semiImplicitSlipDisplacementLaplacianSolidBody`:** Semi-implicit slip fvMotionSolver based on Laplacian smoothing of the point displacement field. Solid-body rotation is added on top of the deformed mesh
  
**1.3. `semiImplicitSlipVelocityLaplacian`:** Semi-implicit slip fvMotionSolver based on Laplacian smoothing of the point velocity field

**1.4. `semiImplicitSlipVelocityLaplacianSolidBody`:** Semi-implicit slip fvMotionSolver based on Laplacian smoothing of the point velocity field. Solid-body rotation is added on top of the deformed mesh


### 2. src/dirichletBC
 
**2.1. `displacementDirichletBC`:** General Dirichlet boundary condition for point displacement field.

**2.2. `velocityDirichletBC`:** General Dirichlet boundary condition for point velocity field.
  
**2.3. `bumpDirichletBC`:** _Ad-hoc_ Dirichlet boundary conditions used in the bump tutorial case.

**2.4. `simpleKaplan`:** _Ad-hoc_ Dirichlet boundary conditions used in the simpleKaplan tutorial case.



### 3. src/rotatingBC
 
**3.1. `displacemenRotatingBC`:** General rotating boundary condition for point displacement field.

**3.2. `velocityRotatingBC`:** General rotating boundary condition for point velocity field.



### 4. src/rotatingSlip
 
**4.1. `fvPatchField`:** Rotating slip boundary condition for fvPatchField (boundary faces).

**4.2. `pointPatchField`:** Rotating slip boundary condition for pointPatchField (boundary points).


## Reference
[1]. S. Salehi, H. Nilsson, An OpenFOAM mesh motion library for CFD of transient operation of Kaplan turbines, _submitted to Computer Physics Communications_.

[2]. S. Salehi, H. Nilsson, E. Lillberg, N. Edh, Development of a novel numerical framework in OpenFOAM to simulate kaplan turbine transients, IOP Conference Series: Earth and Environmental Science 774 (1) (2021) 012058. [doi:10.1088/1755-1315/774/1/012058](https://www.doi.org/10.1088/1755-1315/774/1/012058).
