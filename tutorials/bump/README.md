# bump

## General information

The performance of the proposed semi-implicit slip algorithm is assessed using a two-dimensional mesh deformation test case. 
Figure below illustrates the configuration of this verification case study, which is a simple structured hexahedral mesh on top of 
a bump with a smooth curved surface. The left (green) boundary is moving with a constant speed from left to right while the right (blue) boundary is fixed.
The points are supposed to slip on the upper and lower boundaries (red). Although the designed configuration seems rather simple, 
it is quite challenging for a motion displacement solver due to the slipping points on the curved surface. 
This simple verification case study can be seen as a cheap test case that mimics the complex mesh deformation and slip procedure in Kaplan turbine
transients. The curved lower surface may be considered as the hub or shroud of a Kaplan turbine, on which points should be able to slip smoothly
in a transient sequence.

![test](https://user-images.githubusercontent.com/103576002/169289486-994da310-fe1f-4f96-a589-a02783e3d55d.png)


The total length of the channel is 1.2 m, whereas the bump length is l=1 m. The channel and bump heights are H=0.4 m and h=0.1 m, respectively. The left boundary moves with a constant speed of U=0.1 m/s. The `codedFixedValue` type boundary condition is used for the point displacement field of the left boundary to move the points with a constant horizontal velocity of U=0.1 m/s while accommodating their vertical position to the changes in the boundary length. The right boundary points are fixed using a `fixedValue` condition. The top boundary is a flat surface. Consequently, the `fixedNormalSlip` boundary condition is a proper choice for this patch. However, the lower patch, where the points should be able to slip on the bump is the main challenge. Three different options are examined here, namely, slip (`slip`), slip on a prespecified surface (`surfaceSlipDisplacement`), and the developed semi-implicit slip. The first two options are only boundary conditions for the `pointDisplacement` field, whereas the semi-implicit slip approach is a combination of a new mesh motion solver and a new boundary condition, described in the previous section.

Figure below displays the mesh motion solution for the employed options at different times (every 0.65~s). The slip condition combined with the displacement-based solver is very sensitive to small errors, and tiny distortions in the mesh grow drastically which leads to early destruction of the mesh. The lower boundary wrinkles soon after the start of the simulation and the mesh destroys with negative volumes.

The most appropriate originally available option in OpenFOAM seems to be the slip-on-surface boundary condition (`surfaceSlipDisplacement`). The points are always projected to the specified surface. Thus, they will always remain on the surface during the morphing process and the surface geometry is perfectly preserved. However, as previously described, this is an explicit correction and does not affect the linear system of the mesh motion equations. Although the points remain on the lower surface, the internal points do not follow its curvature and hit the surface which destroys the mesh with negative volumes as a consequence.

In contrast, the developed semi-implicit slip algorithm provides a very stable and smooth mesh deformation. Not only do the mesh points stay on the lower curved surface, but the internal points feel the presence of the bump and follow its curvature. The points on the left-hand side of the bump move upwards, while those on the right-hand side move downwards.

![image](https://user-images.githubusercontent.com/103576002/169289211-f1ff468d-dca0-4230-b439-52c33390b05b.png)


## Runnig and postprocessing the case

The `Allrun` script, in the Simulation folder, runs all the cases. The PostProcessing folder contains a MATLAB script that extracts the mesh quality measures and plots the curves as presented in the paper.
