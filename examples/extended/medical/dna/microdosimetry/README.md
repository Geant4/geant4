\page Examplemicrodosimetry Example microdosimetry

\author S. Incerti (a, *), , H. Tran (a, *), V. Ivantchenko (b), M. Karamitros\n
a. LP2i, IN2P3 / CNRS / Bordeaux 1 University, 33175 Gradignan, France \n
b. G4AI Ltd., UK
* e-mail: incerti@lp2ib.in2p3.fr or tran@lp2ib.in2p3.fr \n

## INTRODUCTION.

The microdosimetry example shows how to use Geant4 and Geant4-DNA physics models
in different regions of the geometry.

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

The geometry is a 10-micron side cube (World) made of liquid water (G4_WATER
material) containing a 2 micron-thick slice (along X) of water (Target).

Particles are shot from the World volume.

The variable density feature of materials is illustrated in DetectorConstruction.
The material can be changed directly in microdosimetry.in macro file.

## SET-UP

Make sure $G4LEDATA points to the low energy electromagnetic data files.

## HOW TO RUN THE EXAMPLE

In interactive mode, run:

```
./microdosimetry
```

In batch, the macro microdosimetry.in can be used. It shows how to shoot different
particle types.

## PHYSICS

The PhysicsList uses Geant4 Physics in the World region and Geant4-DNA Physics
in the Target region.

1) Geant4 Physics in the World is selected via the command:

```
/dna/test/addPhysics X
```

where X is any EM physics list, such as emstandard_opt4 (see PhysicsList.cc).

2) Geant4-DNA activator is used in the regionTarget region using:

```
/process/em/AddDNARegion regionTarget DNA_OptY
```

where Y = 0, 2, 4, or 6.

3) In addition to 1) or 2), to enable radioactive decay, one can use:

```
/dna/test/addPhysics raddecay
```

4) Warning regarding ions: when the incident particle type is ion
(/gun/particle ion), specified with Z and A numbers (/gun/ion A Z),
the Rudd ionisation extended model is used. The particles are tracked
by default down to 0.5 MeV/u. This tracking cut can be bypassed using :

```
/dna/test/addIonsTrackingCut false
```

## SIMULATION OUTPUT AND RESULT ANALYSIS

The output results consists in a dna.root file, containing for each simulation step:
- the type of particle for the current step
- the type of process for the current step
- the step PostStepPoint coordinates (in nm)
- the energy deposit along the current step (in eV)
- the step length (in nm)
- the total energy loss along the current step (in eV)
- the kinetic energy at PreStepPoint (in eV)
- the cos of the scattering angle
- the event ID
- the track ID
- the parent track ID
- the step number

This information is extracted from the SteppingAction class.

The ROOT file can be easily analyzed using for example the provided ROOT macro
file plot.C; to do so :
* be sure to have ROOT installed on your machine
* be sure to be in the directory containing the ROOT files created by microdosimetry
* copy plot.C into this directory
* from there, launch ROOT by typing root
* under your ROOT session, type in : .X plot.C to execute the macro file
* alternatively you can type directly under your session : root plot.C

The naming scheme on the displayed ROOT plots is as follows (see SteppingAction.cc):

- particles: \n
gamma: 0 \n
e-: 1 \n
proton: 2 \n
hydrogen: 3 \n
alpha: 4 \n
alpha+: 5 \n
helium: 6 \n
\n
- processes: \n

Capture: 1 \n
(only if one uses G4EmmicrodosimetryActivator in PhysicsList)

e-_G4DNAElectronSolvation: 10 \n
e-_G4DNAElastic: 11 \n
e-_G4DNAExcitation: 12 \n
e-_G4DNAIonisation: 13 \n
e-_G4DNAAttachment: 14 \n
e-_G4DNAVibExcitation: 15 \n
msc: 110 \n
CoulombScat: 120 \n
eIoni: 130 \n \n

proton_G4DNAElastic: 21 \n
proton_G4DNAExcitation: 22 \n
proton_G4DNAIonisation: 23 \n
proton_G4DNAChargeDecrease: 24 \n
msc: 210 \n
CoulombScat: 220 \n
hIoni: 230 \n
nuclearStopping: 240 \n \n

hydrogen_G4DNAElastic: 31 \n
hydrogen_G4DNAExcitation: 32 \n
hydrogen_G4DNAIonisation: 33 \n
hydrogen_G4DNAChargeIncrease: 35 \n \n

alpha_G4DNAElastic: 41 \n
alpha_G4DNAExcitation: 42 \n
alpha_G4DNAIonisation: 43 \n
alpha_G4DNAChargeDecrease:	44 \n
msc: 410 \n
CoulombScat: 420 \n
ionIoni: 430 \n
nuclearStopping: 440 \n \n

alpha+_G4DNAElastic: 51 \n
alpha+_G4DNAExcitation:	52 \n
alpha+_G4DNAIonisation: 53 \n
alpha+_G4DNAChargeDecrease: 54 \n
alpha+_G4DNAChargeIncrease: 55 \n
msc: 510 \n
CoulombScat: 520 \n
hIoni: 530 \n
nuclearStopping: 540 \n

helium_G4DNAElastic: 61 \n
helium_G4DNAExcitation: 62 \n
helium_G4DNAIonisation: 63 \n
helium_G4DNAChargeIncrease: 65 \n \n

GenericIon_G4DNAIonisation: 73 \n
msc: 710 \n
CoulombScat: 720 \n
ionIoni: 730 \n
nuclearStopping: 740 \n \n

phot: 81 \n
compt: 82 \n
conv: 83 \n
Rayl: 84 \n

---------------------------------------------------------------------------

Should you have any enquiry, please do not hesitate to contact:
incerti@cenbg.in2p3.fr or tran@lp2ib.in2p3.fr
