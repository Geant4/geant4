\page Examplemfp Example mfp

\author S. Incerti et al. (a, *) \n
a. LP2i, IN2P3 / CNRS / Bordeaux University, 33175 Gradignan, France \n
* e-mail:incerti@lp2ib.in2p3.fr

## INTRODUCTION.

The mfp example shows how to calculate mean free path of particles
in liquid water using the Geant4-DNA physics processes and models.

It has been adapted from the spower and TestEm14 examples.

This example is provided by the Geant4-DNA collaboration.

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

The geometry is a 1 m radius sphere of liquid water (G4_WATER
material). Particles are shot along x from the sphere centre.

Radius of the sphere, physics constructor, primary particle type and
energy can be controlled by the mfp.in macro file.

## SET-UP

Make sure G4LEDATA points to the low energy electromagnetic data files.

The code can be compiled with cmake.

It works in MT mode.

## HOW TO RUN THE EXAMPLE

Use:

```
./mfp mfp.in
```

The mfp.in macro allows a full control of the simulation.

The computation of MFP and other quantities is performed in the
SteppingAction::UserSteppingAction method.

The histo.in macro shows how to display several quantities
(energy spectrum, scattering angle along x) of primary and secondaries.

## PHYSICS

G4EmDNAPhysics* constructors are used.

## SIMULATION OUTPUT AND RESULT ANALYSIS

The accuracy of results may depend on incident statistics as well as
on number of steps specified in the SteppingAction::UserSteppingAction
method.

The output results consist in a text file (mfp.txt), containing:
- energy of incident particles (in eV)
- mfp (in nm)
- rms (i.e. standard deviation) on mfp (in nm)

Otherwise you may use histo.in to generate ROOT histograms of the
other quantities.