     =========================================================
 		  Geant4 - range example
     =========================================================

                                README file
                          ----------------------

                           CORRESPONDING AUTHOR

S. Incerti et al. (a, *)
a. LP2i, IN2P3 / CNRS / Bordeaux University, 33175 Gradignan, France
e-mail: incerti@lp2ib.in2p3.fr

---->0. INTRODUCTION.

The range example shows how to calculate range of electrons
in liquid water using the Geant4-DNA physics processes and models.

It is adapted from svalue.

This example is provided by the Geant4-DNA collaboration.

These processes and models are further described at:
http://geant4-dna.org

Any report or published results obtained using the Geant4-DNA software shall
cite the following Geant4-DNA collaboration publications:
Med. Phys. 51 (2024) 5873–5889
Med. Phys. 45 (2018) e722-e739
Phys. Med. 31 (2015) 861-874
Med. Phys. 37 (2010) 4692-4708
Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178

---->1. GEOMETRY SET-UP.

The geometry is a 1 m radius sphere of liquid water (G4_WATER
material). Particles are shot randomly from the sphere centre.

Radius of the sphere, physics constructor and energy can be
controlled by the range.in macro file.

---->2. SET-UP

Make sure G4LEDATA points to the low energy electromagnetic data files.

The code can be compiled with cmake.

It works in MT mode.

---->3. HOW TO RUN THE EXAMPLE

In interactive mode, run:

./range range.in

The range.in macro allows a full control of the simulation.

The histo.in macro is also provided for the creation of histograms.

The computation of ranges is performed in the
TrackingAction::PostUserTrackingAction method. The computation
for electrons uses the primary particle track length and the computation
for incident particles undergoing Geant4-DNA charge exchange processes, such
as protons, hydrogen, alpha particles and their charge states, is specific.

---->4. PHYSICS

You can select Geant4-DNA physics constructor in range.in.

A tracking cut can be applied if requested.

---->5. SIMULATION OUTPUT AND RESULT ANALYSIS

The output results consist in a text file (range.txt), containing :
- energy of incident particles (in eV)
- track length of primary particle (in nm)
- rms of track length of primary particle (in nm)
- projected length of primary particle (in nm)
- rms of projected length of primary particle (in nm)
- penetration of primary particle (in nm)
- rms of penetration of primary particle (in nm)

Note: rms values correspond to standard deviation.

Results in this file can be displayed using the ROOT macro plot.C:
root plot.C
