\page ExampleTestEm6 Example TestEm6 

  This example is intended to test the processes of gamma conversion
  to a pair of muons and annihilation of positrons with atomic
  electrons to a pair of muons.
	
## GEOMETRY DEFINITION
 
  The geometry consists of a single block of a homogenous material.
     
  Two parameters define the geometry :
     - the material of the box,
     - the (full) size of the box.
  The default is 500 m of iron.
             
  In addition a transverse uniform magnetic field can be applied.
 
  The default geometry is constructed in DetectorConstruction class,
  but all of the above parameters can be changed interactively via
  the commands defined in the DetectorMessenger class.
 	
## PHYSICS LIST

  Physics Lists are based on modular design. Several modules are
  instantiated:
  1. Transportation
  2. EM physics
  3. Decays
  4. StepMax - for step limitation

  The electromagnetic physics is chosen from one of the Geant4 EM
  physics constructors in the physics_list library.

  Cross sections can be enhanced (see below).
 	 
## AN EVENT : THE PRIMARY GENERATOR
 
  The primary kinematic consists of a single particle which hits the
  block perpendicular to the input face. The type of the particle
  and its energy are set in the PrimaryGeneratorAction class, and can
  changed via the G4 build-in commands of G4ParticleGun class (see
  the macros provided with this example).
  The default is a Gamma of 100 TeV.
     
  In addition one can choose randomly the impact point of the incident
  particle. The corresponding interactive command is built in
  PrimaryGeneratorMessenger class.
             
  A RUN is a set of events.
 				
## VISUALIZATION
 
  The Visualization Manager is set in the main().
  The initialisation of the drawing is done via the command
```
> /control/execute vis.mac
```
 	
  The detector has a default view which is a longitudinal view of the box.
 
  The tracks are drawn at the end of event, and erased at the end of run.
  Optionally one can choose to draw all particles, only the charged ones,
  or none. This command is defined in EventActionMessenger class.

## PHYSICS DEMO

  The particle's type and the physics processes which will be available
  in this example are set in PhysicsList class.

  In addition a build-in interactive command (/process/inactivate procname)
  allows to activate/inactivate the processes one by one.

  The threshold for producing secondaries can be changed.
  eg:
```
/run/particle/setCut 100 micrometer
/run/initialize
```

  To visualize the GammaConversionToMuons :
```
/control/execute run01.mac
/control/execute vis.mac
/run/beamOn
```

  To visualize the AnnihiToMuPair :
```
/control/execute run11.mac
/control/execute vis.mac
/run/beamOn
```

  Other macros:
  - run02.mac: the final state of the GammaConversionToMuons
  - run12.mac: test on carbon target with biasing of cross section

## HOW TO START ?
 
  - Execute Test  in 'batch' mode from macro files
```
% ./TestEm6    run01.mac
```
 		
  - Execute Test  in 'interactive mode' with visualization
```
% ./TestEm6 
	....
Idle> type your commands
	....
Idle> exit
```
 
## HOW TO INCREASE STATISTICS ON gamma -> mu+mu- ?
 
  The processes of gamma -> mu+mu-  and e+e- -> mu+mu-
  have a low cross section but can be important
  for leakage through thick absorbers and calorimeters.
  Straight forward simulation will be quite time consuming.
  To make the processes more visible, the cross section can be
  artificially increased by some factor (here 1000)
  using the commands  (only effective after  /run/initialize)

```
/testem/phys/SetGammaToMuPairFac  1000
/testem/phys/SetAnnihiToMuPairFac 1000
```
 
	
## HISTOGRAMS
 
  Testem6 produces 6 histograms, h1 - h6, which illustrate
  the final state of the GammaConversionToMuons. The histograms are produced
  with run02.mac and can be displayed with the ROOT macro plotHisto.C.

  The remaining histograms h7 - h16 show various cross sections and h17 the ratio
  of eeToHadr/eeToMu, see their definitions in RunAction.cc

  By default the histograms are saved as testem6.root

  The format of the histogram file can be : root (default), xml, csv,
  by selecting the analysis manager default file type in RunAction.cc
  