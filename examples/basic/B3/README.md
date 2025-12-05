\page ExampleB3 Example B3

 This example simulates schematically a Positron Emitted Tomography system.

## GEOMETRY DEFINITION

   The support of gamma detection are scintillating crystals. A small number
   of such crystals are optically grouped in a matrix of crystals. In
   this example, individual crystals are not described; only the matrix of
   crystals is and it is still called 'Crystal' hereafter.

   Crystals are circularly arranged to form a ring. Few rings make up the full
   detector (gamma camera). This is done by positionning Crystals in
   Ring with an appropriate rotation matrix. Several copies of Ring are
   then placed in the full detector.

   The head of a patient is schematised as a homogeneous cylinder of brain
   tissue, placed at the center of full detector.

   The Crystal material, Lu2SiO5, is not included in the G4Nist database.
   Therefore, it is explicitly built in DefineMaterials().

## PHYSICS LIST

   The physics list contains standard electromagnetic processes and the
   radioactiveDecay module for GenericIon. It is defined in the B3::PhysicsList
   class as a Geant4 modular physics list with registered physics builders
   provided in Geant4:
   - G4DecayPhysics - defines all particles and their decay processes
   - G4RadioactiveDecayPhysics - defines radioactiveDecay for GenericIon
   - G4EmStandardPhysics - defines all EM standard processes

   This physics list requires data files for:
   - low energy electromagnetic processes which path is defined via
     the G4LEDATA envirnoment variable
   - data files for nuclides properties which path is defined via
     the G4ENSDFSTATEDATA envirnoment variable
   - radioactive decay hadronic processes which path is defined via
     the G4RADIOACTIVEDATA envirnoment variable.

   See more on installation of the datasets in
   <a href="http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides
                                         /InstallationGuide/html/ch03s03.html">
   Geant4 Installation Guide, Chapter 3.3: Note On Geant4 Datasets </a>.

## ACTION INITALIZATION

   B3a::ActionInitialization class (see also B3b::ActionInitialization) instantiates and registers to Geant4 kernel  all user action classes.

   While in sequential mode the action classes are instatiated just once,
   via invoking the method:
      B3a::ActionInitialization::Build()
      (see also B3b::ActionInitialization::Build)
   in multi-threading mode the same method is invoked for each thread worker
   and so all user action classes are defined thread-local.

   A run action class is instantiated both thread-local
   and global that's why its instance is created also in the method
      B3a::ActionInitialization::BuildForMaster()
      (see also B3b::ActionInitialization::Build)
   which is invoked only in multi-threading mode.

   Beta decay of Fluor generates a neutrino. One wishes not to track this
   neutrino; therefore one kills it immediately, before created particles
   are put in a stack. This is done via the G4RunManager::SetDefaultClassification()
   call in the Build() function.

## PRIMARY GENERATOR

   The default particle beam is an ion (F18), at rest, randomly distributed
   within a zone inside a patient and is defined in
   B3::PrimaryGeneratorAction::GeneratePrimaries().
   The type of a primary particle can be changed with G4ParticleGun commands
   (see run2.mac).

## DETECTOR RESPONSE : scorers

   A 'good' event is an event in which an identical energy of 511 keV is
   deposited in two separate Crystals. A count of the number of such events
   corresponds to a measure of the efficiency of the PET system.
   The total dose deposited in a patient during a run is also computed.

   Scorers are defined in B3::DetectorConstruction::ConstructSDandField(). There are
   two G4MultiFunctionalDetector objects: one for the Crystal (EnergyDeposit),
   and one for the Patient (DoseDeposit)

   The scorers hits are saved in form of ntuples in a Root file using Geant4
   analysis tools. This feature is activated in the main() function with instantiating
   G4TScoreNtupleWriter.

   Two variants of accumulation event statistics in a run are demonstrated
   in this example:

   B3a:

   At the end of event, the values acummulated in B3a::EventAction are passed
   in B3a::RunAction and summed over the whole run (see B3a::EventAction::EndOfevent()).
   In multi-threading mode the data accumulated in G4Accumulable objects per
   workers is merged to the master in B3a::RunAction::EndOfRunAction() and the final
   result is printed on the screen.

   G4Accumulable<> type instead of G4double and G4int types is used for the B3a::RunAction
   data members in order to facilitate merging of the values accumulated on workers
   to the master.  Currently the accumulables have to be registered to G4AccumulablesManager
   and G4AccumulablesManager::Merge() has to be called from the users code. This is planned
   to be further simplified with a closer integration of G4Accumulable classes in
   the Geant4 kernel next year.

   B3b:

   B3b::Run::RecordEvent(), called at end of event, collects informations
   event per event from the hits collections, and accumulates statistic for
   B3b::RunAction::EndOfRunAction().
   In addition, results for dose are accumulated in a
   standard floating-point summation and using a new lightweight statistical
   class called G4StatAnalysis. The G4StatAnalysis class records four values:
   (1) the sum, (2) sum^2, (3) number of entries, and (4) the number of entries
   less than mean * machine-epsilon (the machine epsilon is the difference
   between 1.0 and the next value representable by the floating-point type).
   From these 4 values, G4StatAnalysis provides the mean, FOM, relative error,
   standard deviation, variance, coefficient of variation, efficiency, r2int,
   and r2eff.

   In multi-threading mode the statistics accumulated per workers is merged
   to the master in B3b::Run::Merge().

<hr>

The following paragraphs are common to all basic examples

## VISUALISATION

   The visualization manager is set via the G4VisExecutive class
   in the main() function in exampleB3.cc.
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
   in the main() function in exampleB3.cc

   The selection of the user command interface is then done automatically
   according to the Geant4 configuration or it can be done explicitly via
   the third argument of the G4UIExecutive constructor (see exampleB4a.cc).

   The gui.mac macros are provided in examples B2, B4 and B5. This macro
   is automatically executed if Geant4 is built with any GUI session.
   It is also possible to customise the icons menu bar which is
   demonstrated in the icons.mac macro in example B5.

## HOW TO RUN

   - Execute exampleB3a in the 'interactive mode' with visualization
```
% ./exampleB3a
and type in the commands from run1.mac line by line:
Idle> /control/verbose 2
Idle> /tracking/verbose 2
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

   - Execute exampleB3a in the 'batch' mode from macro files
   (without visualization)
```
% ./exampleB3a run2.mac
% ./exampleB3a exampleB3.in > exampleB3.out
```
