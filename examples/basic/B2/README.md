\page ExampleB2 Example B2

This example simulates a simplified fixed target experiment.

## GEOMETRY DEFINITION

   The setup consists of a target followed by six chambers of increasing
   transverse size at defined instances from the target. These chambers are
   located in a region called the Tracker region.
   Their shape are cylinders, constructed as simple cylinders
   (in B2a::DetectorConstruction) and as parametrised volumes
   (in B2b::DetectorConstruction), see also B2b::ChamberParameterisation class.

   In addition, a global, uniform, and transverse magnetic field can be
   applied using G4GlobalMagFieldMessenger, instantiated in
   B2a::DetectorConstruction::ConstructSDandField with a non zero field value,
   or via interactive commands.
   For example:
```
/globalField/setValue 0.2 0 0 tesla
```
   An instance of the B2::TrackerSD class is created and associated with each
   logical chamber volume (in B2a) and with the one G4LogicalVolume associated
   with G4PVParameterised (in B2b).

   One can change the materials of the target and the chambers
   interactively via the commands defined in B2a::DetectorMessenger
   (or B2b::DetectorMessenger). For example:
```
/B2/det/setTargetMaterial G4_WATER
/B2/det/setChamberMaterial G4_Ar
```

## PHYSICS LIST

   The particle's type and the physic processes which will be available
   in this example are set in the FTFP_BERT physics list. This physics list
   requires data files for electromagnetic and hadronic processes.
   See more on installation of the datasets in
   <a href="http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/InstallationGuide/html/ch03s03.html">
   Geant4 Installation Guide, Chapter 3.3: Note On Geant4 Datasets </a>.
   The following datasets: G4LEDATA, G4LEVELGAMMADATA, G4SAIDXSDATA and
   G4ENSDFSTATEDATA are mandatory for this example.

   In addition, the build-in interactive command:
```
/process/(in)activate processName
```
   allows the user to activate/inactivate the processes one by one.

## ACTION INITALIZATION

   A newly introduced class, B2::ActionInitialization,
   instantiates and registers to Geant4 kernel all user action classes.

   While in sequential mode the action classes are instatiated just once,
   via invoking the method:
     B2::ActionInitialization::Build()
   in multi-threading mode the same method is invoked for each thread worker
   and so all user action classes are defined thread-local.

   A run action class is instantiated both thread-local
   and global that's why its instance is created also in the method
      B2::ActionInitialization::BuildForMaster()
   which is invoked only in multi-threading mode.

## PRIMARY GENERATOR

   The primary generator action class employs the G4ParticleGun.
   The primary kinematics consists of a single particle which starts close
   to the world boundary and hits the target perpendicular to the entrance
   face. The type of the particle and its energy can be changed via the
   Geant4 built-in commands of the G4ParticleGun class.

## RUNS and EVENTS

   A run is a set of events.

   The user has control:
      - at Begin and End of each run (class B2::RunAction)
      - at Begin and End of each event (class B2::EventAction)
      - at Begin and End of each track (class TrackingAction, not used here)
      - at End of each step (class SteppingAction, not used here)

   The event number is written to the log file every requested number
   of events in B2::EventAction::BeginOfEventAction() and
   B2::EventAction::EndOfEventAction().
   Moreover, for the first 100 events and every 100 events thereafter
   information about the number of stored trajectories in the event
   is printed as well as the number of hits stored in the G4VHitsCollection.

   The run number is printed at B2::RunAction::BeginOfRunAction(), where the
   G4RunManager is also informed how to SetRandomNumberStore for storing
   initial random number seeds per run or per event.

## USER LIMITS

   This example also illustrates how to introduce tracking constraints
   like maximum step length, minimum kinetic energy etc. via the G4UserLimits
   class and associated G4StepLimiter and G4UserSpecialCuts processes.
   See B2a::DetectorConstruction (or B2b::DetectorConstruction).

   The maximum step limit in the tracker region can be set by the interactive
   command (see B2a::DetectorMessenger, B2b::DetectorMessenger classes).
   For example:

```
/B2/det/stepMax 1.0 mm
```

## DETECTOR RESPONSE

   A HIT is a step per step record of all the information needed to
   simulate and analyse the detector response.

   In this example the Tracker chambers are considered to be the detector.
   Therefore, the chambers are declared 'sensitive detectors' (SD) in
   the B2a::DetectorConstruction (or B2b::DetectorConstruction) class.
   They are associated with an instance of the B2::TrackerSD class.

   Then, a Hit is defined as a set of 4 informations per step, inside
   the chambers, namely:
     - the track identifier (an integer),
     - the chamber number,
     - the total energy deposit in this step, and
     - the position of the energy deposit.

   A given hit is an instance of the class B2::TrackerHit which is created
   during the tracking of a particle, step by step, in the method
   B2::TrackerSD::ProcessHits(). This hit is inserted in a HitsCollection.

   The HitsCollection is printed at the end of each event (via the method
   B2::TrackerSD::EndOfEvent()), under the control of the command:
```
/hits/verbose 2
```

<hr>

 The following paragraphs are common to all basic examples

## VISUALISATION

   The visualization manager is set via the G4VisExecutive class
   in the main() function in exampleB2a.cc (or exampleB2b.cc).
   The initialisation of the drawing is done via a set of /vis/ commands
   in the macro vis.mac. This macro is automatically read from
   the main function when the example is used in interactive running mode.

   By default, vis.mac opens the default viewer (/vis/open).
   This chooses a graphics system (in order of priority):
   - by argument in G4VisExecutive construction.
   - by environment variable, G4VIS_DEFAULT_DRIVER.
   - by information in ~/.g4session.
   - by mode (batch/interactive) and if interactive, by your build flags.

   The user can change the initial viewer
   - with environment variable G4VIS_DEFAULT_DRIVER. The format is
     ```
     <graphics-system> [<window-size-hint>]
     ```
     Set this, e.g:
     - (bash) export G4VIS_DEFAULT_DRIVER=TSG
     - (tcsh) setenv G4VIS_DEFAULT_DRIVER OI
       - The window-size-hint can optionally be added, e.g:
       - (bash) export G4VIS_DEFAULT_DRIVER="RayTracerQt 1000x1000-0+0"
   - on the command line, precede the app invocation, e.g:
     - ```
       G4VIS_DEFAULT_DRIVER=Vtk ./<application-name>
       ```
   - with ~/.g4session.

   For other suggestions for G4VIS_DEFAULT_DRIVER (see list of registered
   graphics systems printed at the start):
   - DAWNFILE: to create a .prim file suitable for viewing in DAWN.
   - VRML2FILE: to create a .wrl file suitable for viewing in a VRML viewer.
   - "TSG_OFFSCREEN 1200x1200": to create an image file with TSG.
     - See the tsg_offscreen.mac in examples/basic/B5 for more commands
       to change the file format, file name, picture size, etc.

   See "Choosing a graphics viewer" in the Application Guide for details.

   Of course you can change the viewer by editing the /vis/open line in vis.mac.

   Also, after the initial viewer opens, you may open a different viewer by typing
   on the command line, e.g:
```
/vis/open DAWNFILE
```
   or
```
/vis/open RayTraceQt
```
   (if you are using the Qt GUI).

   The view parameters of the existing viewer are copied.

   The DAWNFILE and similar drivers are always available
   (since they require no external libraries), but the OGL driver requires
   that the Geant4 libraries have been built with the OpenGL option.

## USER INTERFACES

   The user command interface is set via the G4UIExecutive class
   in the main() function in exampleB2a.cc

   The selection of the user command interface is then done automatically
   according to the Geant4 configuration or it can be done explicitly via
   the third argument of the G4UIExecutive constructor (see exampleB2a.cc).

   The gui.mac macros are provided in examples B2, B4 and B5. This macro
   is automatically executed if Geant4 is built with any GUI session.
   It is also possible to customise the icons menu bar which is
   demonstrated in the icons.mac macro in example B5.

## HOW TO RUN

   - Execute exampleB2a in the 'interactive mode' with visualization
```
% exampleB2a
and type in the commands from run1.mac line by line:
Idle> /tracking/verbose 1
Idle> /run/beamOn 1
Idle> ...
Idle> exit
```
      or
```
Idle> /control/execute run1.mac
....
Idle> exit
```

   - Execute exampleB2a in the 'batch' mode from macro files
   (without visualization)
```
% exampleB2a run2.mac
% exampleB2a exampleB2.in > exampleB2.out
```
