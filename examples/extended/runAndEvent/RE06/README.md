\page ExampleRE06 Example RE06

 This example simulates three simplified sandwitch calorimeters.
 The main features demonstrated in this example are :

-# Utilizing a concrete run class derived from G4Run base class for
   accumulating physics quantities for a run
-# Changing calorimeter geometries without re-building a world volume
-# Defining geometrical regions and setting production thresholds
   for each region
-# Demonstrating the use of primitive scorer and filter classes without
   implementing sensitive detector class
-# Demonstrating the use of parallel scoring geometry and associating
   parallel world scoring process
-# Measuring the timing spent for each region, both for all particle
   types and for e+/e-

 It was moved in extended examples from novice/N07 with removal of
 novice examples. 

 <i> Note: Since this example utilizes its own RE06SteppingVerbose for the
 timing measurement, the user cannot get the ordinary verbosity with
 /tracking/verbose. </i>

## Utilizing a concrete run class derived from G4Run base class for accumulating physics quantities for a run

 G4Run is a class the user can inherit and create his/her own concrete
 class for accumulating information useful to him/her. It has a virtual
 method RecordEvent(const G4Event*), which will be invoked by G4RunManager
 at the end of processing each event. By implemeting this method in the
 user'r concrete run class, he/she can store information associating with
 G4Event class itself and hits collections attached with G4Event. In this
 example, RE06Run is the class derived from G4Run. In the method
 RE06Run::RecordEvent(const G4Event*), in addition to counting the
 number of events, all hits collections are accessed to accumulate
 energy depositions, step lengths and number of steps.

 In case the user create his/her own run class, an object of this class
 must be instantiated in the method GenerateRun() of his/her concrete
 class derived from G4UserRunAction base class. The pointer to this run
 object must be returned by this method. In this example, RE06RunAction
 is the class which instantiating RE06Run class object.  In
 RE06RunAction::EndOfRunAction(const G4Run*) method, RE06Run object
 is analized to output the run summary.
         
## Changing calorimeter geometries without re-building a world volume

 In RE06DetectorConstruction, all solids, logical and physical volumes
 are constructed only once at the first invocation of RE06DetectorConstruction::Constuct() method.
 Positions and number of slices are changed not by re-constructing another
 objects but by modifying data members of already existing objects as
 it is implemented in RE06DetectorConstruction::SetNumberOfLayers(G4int)
 for changing the number of parameterized volumes, and also
 RE06DetectorConstruction::SetSerialGeometry(G4bool) for changing the
 position of placed volumes.
 
## Defining geometrical regions and setting production thresholds for each region

  Setting production thresholds (so-called production cuts) to individual
  region of a detector geometry is the new feature provided by Geant4 5.1
  release. This feature is also called as "Cuts per region".

  Please note that this new feature is supporsed to be used only by the
  users,
   a) who is simulating most complex geometry such as an LHC detector,
   b) and who has enough experience of simulating EM showers in matter.
  We strongly recommend to compare the simulated results of this new
  feature with the results of the same geometry but having uniform 
  production thresholds. Setting completely different cut values for
  individual region may break the coherent and comprehensive accuracy
  of the simulation. Thus such cut values should be carefully optimized
  by the user with comparison with results of uniform cuts.

  In RE06DetectorConstruction::Construct(), Three objects of G4Region
  class are instantiated and set to the logical volumes of each of three
  calorimeter modules. Also, these individual logical volumes are
  registered as "root logical volume" so that all daghter volumes in
  these logical volumes are also affected by the corresponding regions.

## Demonstrating the use of primitive scorer and filter classes without implementing sensitive detector class

 In RE06DetectorConstruction::SetupDetectors() method, concrete classes
 G4PSEnergyDeposit, G4PSNofSecondary, G4PSTrackLength, G4PSNofStep and
 G4PSMinKinEAtGeneration, all of thich are derivalable of G4VPrimitiveScorer,
 are used to define the sensitivity of the calorimeter. All of them are
 registered to G4MultiFunctionalDetector and this detector object is set
 to the logical volume. G4SDParticleFilter is used to define the particle
 type(s) to be scored.

 In RE06Run::RecordEvent() method, the way of retreiving G4THitsMap
 from each primitive scorer via G4HCofThisEvent is demonstrated.
 In RE06RunAction::EndOfRunAction(), Run is summarized with data kept
 in RE06Run class object.

## Demonstrating the use of parallel scoring geometry and associating parallel world scoring process

 G4ParallelWorldScoringProcess is
 assigned to all the particle types. This process invokes sensitive detectors
 (and scorers) defined in the parallel world "ParallelScoringWorld", the
 name of the parallel world which is defined in main() as
 an argument of RE06ParallelWorld constructor.
 
 As implemented in RE06ParallelWorld::SetupGeometry(), the world volume of
 the parallel world is obtained by GetWorld() method as a clone copy of the
 world volume of the mass geometry. The user should not create the world volume.

 RE06ParallelWorld defines three cylindrical volumes, each of them is
 located at the same position as three sandwitch calorimeters defined
 in the mass geometry (RE06DetectorConstruction). Each cylinder is replicated
 in Rho to define 20 layers, and scores the same quantities as the mass geometry.
 These three cylinders are relocated accordingly when the mass geometry is
 modified by RE06DetectorConstruction::SetSerialGeometry().

## Measuring the timing spent for each region, both for all particle types and for e+/e-

 RE06SteppingVerbose class has two G4SliceTimer class objects for each
 detector region. One G4SliceTimer is measuring the time spent by a step
 in a region for all types of particles, and another is measuring for
 e+/e- only. 

 RE06SteppingVerbose::InitializeTimers() is invoked by RE06RunAction::
 BeginOfRunAction(), and checks the number of regions appear in the
 geometry and instantiates the necessary number of timers. Thus, this
 RE06SteppingVerbose class can be used for any kind of geometry the user
 defines without any modification. Given G4VSteppingVerbose is not invoked
 if the verbosity of G4SteppingManager is 0, this verbosity is set to 1.

 NewStep() and StepInfo() are the methods defined in G4VSteppingVerbose
 base class, and they are invoked at the beginning and the end of every
 step, respectively, from G4SteppingManager. Thus, these methods are
 utilized in RE06SteppingVerbose to start/resume and pause the timer.

 RE06SteppingVerbose::Report() method is used by RE06RunAction::
 EndOfRunAction() to get the timing measured.

## Macro files

 - exampleRE06.in \n
    To be used for batch mode. The reference output file is made by this
    macro file.

 - sample.mac \n
    To be used for interactive mode. Issue "/control/execute sample.mac"
    when "Idle>" prompt appears.

 - vis.mac \n
    Setting visualization parameters. This macro file will be called 
    automatically when interactive execution starts.

## UI commands defined in this example
 
- Select Material of the Absorber, Gap, add materials:
```
/RE06/setAbsMat matName
/RE06/setGapMat matName
/RE06/AddMaterial
```
- Set number of layers:
```
/RE06/numberOfLayers nofLayers
```
- Select calorimeters to be placed in serial or parallel:
```
/RE06/serialGeometry True|False
```

See the complete guidance with `/control/manual RE06`.
