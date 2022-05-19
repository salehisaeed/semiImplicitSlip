# kaplanFoam:

## General information
**kaplanFoam** is a mesh motion library in OpenFOAM for CFD simulation of Kaplan turbines transient operations in which the operating condition is changed by simultaneously adjusting the guide vane and runner blade angles. CFD simulations of such transients are highly complex as it involves mesh deformation of the guide vane passage and simultaneous mesh deformation and rigid-body rotation of the runner blade passage. In particular, the mesh motion process requires a proper slipping of mesh points on highly curved surfaces. The current explicit slip boundary conditions in OpenFOAM fail to make large deformations of the mesh without inverting some of the cells. Thus, a robust semi-implicit slip algorithm is developed, based on the Laplacian smoothing methodology, to tackle this issue. The algorithm, which is implemented as a dynamic library in OpenFOAM, includes different mesh motion solvers and boundary conditions. The library is verified using simple, yet relevant, test cases, including a simple Kaplan turbine.

<!-- All the details about the codes and test cases are presented in a submitted manuscipt (to be published). -->

The library is tested and verified with the recent versions of [OpenFOAM](https://www.openfoam.com/), such as, v1912, v2006, v2012, v2106, v2112. Please use the [Github](https://github.com/salehisaeed/kaplanFoam) version to receive the most recent updates .



## Installation

The library is developed based on OpenFOAM ([openfoam.com](https://www.openfoam.com/)). Therefore, a complete installation of OpenFOAM (prefereably a recent version, such as OpenFOAM-v2112, is required.

```bash
pip install foobar
```