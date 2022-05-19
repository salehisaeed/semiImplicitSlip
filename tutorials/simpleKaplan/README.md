# simpleKaplan test case

## General information
The U9 Kaplan turbine model that is chosen as a verification case study in the paper, is a commercial turbine model owned by [Vattenfall](https://www.vattenfall.se/english/) company and thus the geometry cannot be presented as an open-source test case here. Therefore, a simple Kaplan turbine is designed and presented as a verification case study.

One should note that the _simpleKaplan_ case is merely produced to show the performance of the developed library and the geometry is not carefully designed or optimized. The NACA0012 airfoil is chosen as the guide vane profile, while the runner blade profile is a simple thin ellipse with a 20 deg twist in the spanwise direction. The turbine model consists of 24 guide vanes and 7 runner blades. The diameter of the turbine is 0.5 m. The guide vane opening angle is 30 deg.

In order to reduce the computational cost of the verification case study, only one passage of the guide vane and runner regions is considered. The same approach can be simply extended to a full turbine. Additionally, only mesh motion calculations are performed and the flow is not solved. The figure below presents a 3D view of the model assembly.
![domain](https://user-images.githubusercontent.com/103576002/169298377-b88e322d-7f1d-49ab-a92f-138d5eb29d5b.jpeg)

One case see that a part of the shroud surface (the project of the blade tip onto the shroud) is defined as a separate boundary condition. A rotating boundary condition is chosen for that surface. This is not essential for the simulation but it helps the robustness of the dynamic mesh procedure, especially in a real Kaplan turbine where the tip clearance gap is reduced with turbine load reduction.

The time step and write time of the case is designed in such a way that altough the runner is rotating with full speed the results are seen in a stroboscopic view so that the runner blade rotation is clearly visible.

## Running the case
The `Allrun` performs all the necessary preprocessing procedures and runs the case.
