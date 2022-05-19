# kaplanFoam

A mesh motion library in OpenFOAM for CFD simulation of Kaplan turbines transient operation 

During transient operation of Kaplan turbines, the operating condition is changed by simultaneously adjusting the guide vane and runner blade angles. CFD simulations of such transients are highly complex as it involves mesh deformation of the guide vane passage and simultaneous mesh deformation and rigid-body rotation of the runner blade passage. In particular, the mesh motion process requires a proper slipping of mesh points on highly curved surfaces. The current explicit slip boundary conditions in OpenFOAM fail to make large deformations of the mesh without inverting some of the cells. Thus, a robust semi-implicit slip algorithm is in the present work developed, based on the Laplacian smoothing methodology, to tackle this issue. The algorithm, which is implemented as a dynamic library in OpenFOAM, includes different mesh motion solvers and boundary conditions. The library is verified using simple, yet relevant, test cases, including a simple Kaplan turbine.


All the details about the codes and test cases are presented in the submitted manuscipt.

Please use the Github version to get the most recent updates (https://github.com/salehisaeed/kaplanFoam).
