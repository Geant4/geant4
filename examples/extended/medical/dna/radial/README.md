\page Exampleradial Example radial

\author S. Incerti (a, *)\n
a. LP2i, IN2P3 / CNRS / Bordeaux 1 University, 33175 Gradignan, France \n
* e-mail: incerti@lp2ib.in2p3.fr\n

## INTRODUCTION

The radial example shows how to simulate radial dose profiles in liquid water
from incident ions using the Geant4-DNA physics processes and models.

The Geant4-DNA processes and models are further described at:
http://geant4-dna.org

Any report or published results obtained using the Geant4-DNA software shall
cite the following Geant4-DNA collaboration publications: \n
Med. Phys. 51 (2024) 5873–5889 \n
Med. Phys. 45 (2018) e722-e739 \n
Phys. Med. 31 (2015) 861-874   \n
Med. Phys. 37 (2010) 4692-4708 \n
Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178

## GEOMETRY SET-UP

The geometry is a set of cylindrical shells (hollow cylinders) made of liquid water
(G4_WATER material) and aligned along the Z-axis. 

The World is a cylinder wrapping these shells and having the same maximum radius. 
Its material is water and its length can be made larger 
than the length of the cylinders using an offset.

a) The maximum radius of the set can be specified using the following UI command:

```
/radial/setWorldRadius value unit
```

b) The Z length of the shells can be set using:

```
/radial/setCylinderLength value unit
```

c) An offset can be added along Z in front and after the shells using the command:

```
/radial/setWorldOffsetLength value unit
```

d) The thickness of each shell can be specified using the command:

```
/radial/setCylinderThickness value unit
```

Particles are shot from the center of the entrance surface of the offset space, 
placed just before the cylinders along the Z axis. The Z position is negative 
and is equal to the half-length of the World augmented with the offset.

## DATA

Make sure $G4LEDATA points to the low-energy electromagnetic data files.

## HOW TO RUN THE EXAMPLE

In interactive mode, run:

```
./radial
```

In batch, the macro radial.in can be used. It shows how to shoot different
types of ions and how to use Geant4-DNA Physics constructors.

## PHYSICS

The PhysicsList uses Geant4-DNA Physics constructors.

Geant4-DNA Physics constructors can be selected using the command:

```
/radial/addPhysics DNA_OptX
```

where X is 0 to 8 (2, 4 or 6 are recommended).

Comments regarding ions:

- only the ionisation process is considered, and the radial absorbed dose is obtained
considering all energy losses by secondary electrons

- when the incident particle type is ion (/gun/particle ion), specified with Z
and A numbers (/gun/ion A Z), the Rudd ionisation extended model is used.
The particles are tracked by default down to 0.5 MeV/u and undergo below a capture
process. This tracking cut can be bypassed using:

```
/radial/addIonsTrackingCut false
```

## SIMULATION OUTPUT

The output results consist in a ROOT radial.root file, containing an ntuple, named
"radial" and containing the absorbed dose for each cylinder, identified by its
inner radius, per incident ion.

The ROOT macro file plot.C can be used to draw the radial dose profile.
