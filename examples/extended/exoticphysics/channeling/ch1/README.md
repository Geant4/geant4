\page Examplech1 Example ch1

\author Alexei Sytov - INFN Ferrara Division (Italy) \n
 sytov@fe.infn.it

## INTRODUCTION
Example ch1 is an easy demonstration of the minimum requirements necessary
to integrate the G4ChannelingFastSimModel and the G4BaierKatkov model into a project 
in order to simulate the physics of channeling and 
channeling radiation/coherent bremsstrahlung.

This example serves as a guideline for users on how to add this physics 
to their existing Geant4 projects. It includes the minimum necessary options 
to incorporate this physics. Specifically, it requires registering 
G4FastSimulationPhysics in the main routine and 
adding a few lines of code in DetectorConstruction.

## DESCRIPTION
The example is based on the following experiments on channeling [1] and 
channeling radiation [2] in a bent crystal, carried out at Mainz Mikrotron MAMI with 
855 MeV electrons. The experimental validation of G4ChannelingFastSimModel is 
described in [3].

This example includes a bent crystal and a detector positioned behind it.
The incoming beam is set up in macro run.mac.

The example does not include any input of the model or geometry parameters 
from the macro to keep it as straightforward as possible. The output is recorded
into the file results.root. It consists of 
the charged particle distribution at the detector in the x-plane 
(the plane of crystal bending and perpendicular to the crystal planes) as well as
the spectrum of photons arriving to the detector. To build these plots, one has to
open this file in root and use 
```cpp
x_out->Draw() 
```

and 

```cpp
Spectrum->Draw()
```

for the coordinates and the spectrum, respectively.

## REFERENCES

-# A. Mazzolari et al. <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.112.135503">Phys. Rev. Lett. 112, 135503 (2014).</a>
-# L. Bandiera et al. <a href="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.025504">Phys. Rev. Lett. 115, 025504 (2015).</a>
-# A. Sytov et al. <a href="https://link.springer.com/article/10.1007/s40042-023-00834-6"> Journal of the Korean Physical Society 83, 132–139 (2023).</a>
