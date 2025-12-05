\page Examplephasespace Example phasespace

\author S. Incerti (a, *)\n
a. LP2i, IN2P3 / CNRS / Bordeaux 1 University, 33175 Gradignan, France \n
* e-mail: incerti@lp2ib.in2p3.fr \n

## INTRODUCTION.

The phasespace example shows how to simulate particle tracks in a cube of 
liquid water containing a scoring sphere, placed in its centre. The
example creates two phase space files 
- one with the ROOT (psp.root) format,
- one with the GRAS format (GRAS.csv), 
both scoring particles entering the sphere.

## GEOMETRY SET-UP

The geometry is a 100-micron side cube (World) made of liquid water (G4_WATER
material), containing a sphere of 10-micron radius, placed in its center. 

Particles are shot along the Z axis (positive direction), starting at a 
Z = -12 microns from the center of the World.

The World and sphere sizes and material can be changed directly 
in the phasespace.in macro file.

The shooting is controled using GPS commands in phasespace.in.

The variable density feature of materials is illustrated in
DetectorConstruction. The material density can be changed directly in the
phasespace.in macro file.

## SET-UP

Make sure $G4LEDATA points to the low energy electromagnetic data files.

## HOW TO RUN THE EXAMPLE

In interactive mode, run:

```
./phasespace
```

In batch, the macro phasespace.in can be used. It shows how to define dimensions 
and shoot different particle types.

## PHYSICS

The list of available Reference Physics Lists is shown when simulation starts.

All Physics settings are done in the PhysicsList class.
QGSP_BIC_HP is activated by default.

At run time, tracks are killed inside the sphere, 
using in phasespace.in the command :

```
/step/killInsideScorer 1
```

## SIMULATION OUTPUT AND RESULT ANALYSIS

The output results consists in a psp.root file, containing a ntuple, named
"scorer". 

For each particle entering the sphere, when the PreStepPoint is located at 
the geometry boundary, it records :
- the PDG encoding of the particle
- the 3 PreStepPoint coordinates (in microns)
- the 3 PreStepPoint momentum cosine (no unit)
- the kinetic energy at PreStepPoint (in keV)

This information is extracted from the SteppingAction class.

The ROOT file can be easily analyzed using the provided ROOT macro file 
plot.C; to do so :
* be sure to have ROOT installed on your machine
* be sure to be in the directory containing the ROOT file created by phasespace
* from there, launch ROOT by typing root
* under your ROOT session, type in : .X plot.C to execute the macro file
* alternatively you can type directly under your session : root plot.C

The ROOT macro also writes a GRAS.csv file, following the format of GRAS (1). 
Such a file can for ex. be read with the moleculardna Geant4 advanced example (2).

- (1) https://essr.esa.int/project/gras-geant4-radiation-analysis-for-space
- (2) https://moleculardna.org 
