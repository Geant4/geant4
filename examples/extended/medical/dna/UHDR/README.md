\page ExampleUHDR Example UHDR

  This example is provided by the Geant4-DNA collaboration
  (http://geant4-dna.org).

  Any report or published results obtained using the Geant4-DNA software
  shall cite the following Geant4-DNA collaboration publications:\n
  Med. Phys. 45 (2018) e722-e739\n
  Phys. Med. 31 (2015) 861-874\n
  Med. Phys. 37 (2010) 4692-4708\n
  Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178\n


## INTRODUCTION

 This example shows how to activate the mesoscopic model in chemistry and
 combine with IRT-syn model (https://arxiv.org/abs/2409.11993).
 It allows to simulate chemical reactions longtime (beyond 1 us) of post-irradiation.
 Important : This example is provided as a prototype. Some mistakes may be present in the code.
 We would be happy to hear from you — please don’t hesitate to share your feedback and 
 let us know about any difficulties you may encounter.

 To run the example:
```
mkdir UHDR-build
cd UHDR-build
cmake ../pathToExamples/UHDR
make
```

To visualize (only for physical stage)
```
./UHDR
```

In batch mode, the macro beam.in can be used as follows:
```
./UHDR beam.in
```
 or
```
./UHDR beam.in 123
   # 123 is the user's seed number
```

## GEOMETRY DEFINITION


 The world volume is a simple water box 3.2 x 3.2 x 3.2 um3 for 0.01 Gy of cut-off
 absorbed dose and 1.6 x 1.6 x 1.6 um3 for 1 Gy. This example is limited to these geometries.
 The choice of simulation volume is a compromise between a sufficient number of chemical species a
 nd an achievable computation time.

 Two parameters define the geometry :
 - the material of the box for the physical stage is water.
 - for the chemistry stage, the concentration of scavengers in [mole/l]
   is added. This concentration is supposed to have no effect on the
   physical stage. pH is defined as scavengers of H3O^1, OH^-1.
   In this example, we consider that chemical molecules diffuse and react in a
   bounded volume (that is, limited by geometrical boundaries) which is also
   the irradiated water box volume of the physical stage.
   The bouncing of chemical molecules on the volume border is applied
   for both SBS and mesoscopic models.
   The bouncing is not applied for the physical stage.

## PHYSICS LIST

 PhysicsList is Geant4 modular physics list using G4EmDNAPhysics_option2
 and EmDNAChemistry constructors (the chemistry constructor uses the
 Step-by-step method).

## Beamline Dose Cutoff

 A threshold dose value along a beamline in a simulation. 
 This is the total dose per beamline (across all pulses). 
 The concept of the dose cutoff is that once the recorded dose exceeds this threshold, 
 the particle beam is automatically stopped. 

```
 /scorer/Dose/cutoff 0.1 Gy
```

## CHEMISTRY WORLD

 This object is controlled by DetectorContruction. It defines the chemistry volume,
 scavengers and pH of water.

 This fearture can be set by the following commands:
 # pH and Scavenger
```
 /UHDR/env/pH 5.5
```
 # air concentration
```
 /UHDR/env/scavenger O2 21 %
 /UHDR/env/scavenger CO2 0.041 %
 /UHDR/env/scavenger HCO3m 2.4 uM
```
## AN EVENT: PRIMARY GENERATOR

This example utilizes the G4SingleParticleSource.
Each event consists of multiple incident particles.
A large number has been chosen to ensure that the stack remains non-empty until the desired
energy deposition is achieved (which is then converted to a cutoff dose).
With each /run/beamOn command, a group of particles is emitted. The cutoff dose
(dose threshold) determined by users.
The actual dose is calculated based on the real energy deposited in the volume.

## G-value Scorer

There is one G4MultiFunctionalDetector object which computes the
energy deposition and the number of species along time in order to
extract the G-value:
(Number of species X) / (100 eV of deposited energy).

These two macro commands can be used to control the scoring time:
```
/scorer/Gvalues/addTimeToRecord 1 ps
   # user can select time bin to score G values.
/scorer/Gvalues/nOfTimeBins
   # or user can automatically select time bin logarithmically.
```

## PULSE ACTION and INTERPULSE ACTION

The time structure can be activated by implementing a delayed time, Δt,
which is sampled from a beamline raw signal of a measured pulse.
Each delayed time is associated with a primary particle and propagates to its corresponding
primary chemical species induced by this primary particle.
These primary chemical species remain inactive until the virtual simulation time matches
their respective delayed time. This process creates a duration for the primary particle train
where their primary chemical species are activated randomly through an experimental beam
current transformation, named "pulse duration".

This fearture can be set by the following commands:

# time structure
```
/UHDR/pulse/pulseOn true // active the time structure
```
# push structure file
```
/UHDR/pulse/pulseFile 1.4us // push structure file
```
# pulse structure
```
/UHDR/pulse/multiPulse true // active the multi pulse
/UHDR/pulse/pulsePeriod 10 ms // time between two pulses (DIT)
/UHDR/pulse/numberOfPulse 2 // number of pulses
```
 **Important: the multi-pulse function may produce incorrect results with long pulse structures.**
## OUTPUT

G-value

## RELEVANT MACRO COMMANDS AND MACRO FILE

The user macro files are:
- beam.in (default),
- CONV.in (Conventional)
- UHDR.in (Ultra High Dose Rate)
- initialize.in (initialize geo and phys)
- scavengers.in (pH and scavengers are defined)

## REACTION BUILDER

 Reaction lists are collected by builders for specific applications.
 - ChemNO2_NO3ScavengerBuilder is to build the reaction list with NO2-/NO3-.
 - ChemPureWaterBuilder is to build the reaction list with pH.
 - ChemOxygenWaterBuilder is to build the reaction list with ROS.
 - ChemFrickeReactionBuilder is to build the reaction list of Fricke Dosimeter.

## PLOT

The information about all the molecular species is scored in a ROOT
(https://root.cern) ntuple file Dose_xxx.root (xxx is seed number).
The ROOT program plot_time.C
can be used to plot the G values vs time for each species.

Execute plot_time as:
```
root plot_time.C
```

or print G values to scorer.txt
```
root plot_time.C > scorer.txt
```

The results show the molecular species (G values) as a function of
time (ns).

## Periodic Boundary Condition (PBC)

The Periodic Boundary Condition is implemented based on https://github.com/amentumspace/g4pbc
to calculate microdosimetry. The periodic boundary condition (PBC) is used to simulate the
behavior of secondary electrons during the physical stage.
When an electron exits an edge of a cubic volume, it re-enters from the opposite edge.
The PBC helps reduce the edge effects in dose calculations for micrometer-sized volumes


The PBC requires a maximum dose (xxx) to abort the event. This to avoid the high energy of
secondary electrons deposit a large energy inside the micro volume.

```
/scorer/Dose/abortedDose xxx Gy
```

Use the following command to activate or deactivate PBC.

```
/UHDR/Detector/PBC true
```

Contact: H. Tran (tran@lp2ib.in2p3.fr)
CNRS, lp2i, UMR 5797, Université de Bordeaux, F-33170 Gradignan, France