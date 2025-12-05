\page ExampleRE05 Example RE05

 Example RE05 has a simplified collider detector geometry. This example
demonstrates the following features. \n
 It was moved in extended examples from novice/N04 with removal of
novice examples. 

## PYTHIA primary events

 RE05PrimaryGeneratorAction has G4HEPEvtInterface as the generator.
G4HEPEvtInterface accesses to "pythia_event.data", which contains three
events of Higgs generation produced by PYTHIA. "pythia_main.f" is an
example FORTRAN code of PYTHIA for generating this event sample.

## Readout geometry

 RE05DetectorConstruction defines a simplified collider detecor
geometry, a tracker made of cylindrical tubes, a calorimeter made of
cylindrical tubes, and muon trackers made of planes.

 The cylindrical calorimeter is made of tubes of lead and a scintillator.
Energy deposition in the scintillator is accumulated by RE05CalorimeterSD 
sensitive detector, which is assigned to a dedicated parallel world,
RE05CalorimeterParallelWorld, which defines the phi-z cell.

## Physics processes

 The example uses the QBBC physics list, which includes electromagnetic 
and hadronic interactions. 

## Event filtering by the stacking mechanism

 Higgs events in "pythia_event.data" have two lepton pairs produced
by the Higgs decay via Z0. At the first stage of each event, only the
primary muons are tracked without tracking secondaries. then the number
of hits on the muon trackers are examined. At the next stage, only
the primary charged particles are tracked only inside the barrel
tracking area and the isolation of the primary muons are examined.
At the third stage, all particles in the RoI (Region of Interest) along
the isolated muons are tracked. All these examinations are applied in
RE05StackingAction.
  	
## Input macro files

- exampleRE05.in: 

  Read "pythia_event.data" and run 10 events. This macro file
  demonstrates the feature described in the previous section.

  N.B. For release 11.3, this macro is temporarily altered to shoot
  muons by particle gun. Once HepMC3 interface becomes ready, this
  example will be overhauled and this macro file will be updated
  accordingly.

- exampleRE05.EMtest.in and exampleRE05.EMtest.large_N.in

  Alternative macros with shooting muons for checking EM physics.

## How to start
 
- Execute RE05 in 'batch' mode from macro files
```
% ./exampleRE05   exampleRE05.in
```
 		
- Execute RE05 in 'interactive mode' with visualization
```
% ./exampleRE05
....
Idle> type your commands. For instance:
Idle> /run/beamOn 3
....
Idle> exit
```
