\page ExampleB1 Example B1

 This example demonstrates a very simple application where an energy
 deposit is accounted in user actions and their associated objects
 and a dose in a selected volume is calculated.

## GEOMETRY DEFINITION

   The geometry is constructed in the B1::DetectorConstruction class.
   The setup consists of a an envelope of box shape containing two
   volumes: a spherical cone and a trapezoid.

   In this example we use  some common materials materials for medical
   applications. The envelope is made of water and the two inner volumes
   are made from tissue and bone materials.
   The materials are created with the help of the G4NistManager class,
   which allows to build a material from the NIST database using their
   names. Available materials and their compositions can be found in
   <a href="http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides
                                    /ForApplicationDeveloper/html/apas10.html">
   the Geant4 User's Guide for Application Developers, Appendix 10:
   Geant4 Materials Database
   </a>.

## PHYSICS LIST

   The particle's type and the physic processes which will be available
   in this example are set in the QBBC physics list. This physics list
   requires data files for electromagnetic and hadronic processes.
   See more on installation of the datasets in
   <a href="http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides
                                         /InstallationGuide/html/ch03s03.html">
   Geant4 Installation Guide, Chapter 3.3: Note On Geant4 Datasets </a>.
   The following datasets: G4LEDATA, G4LEVELGAMMADATA, G4NEUTRONXSDATA,
   G4SAIDXSDATA and G4ENSDFSTATEDATA are mandatory for this example.

   In addition the build-in interactive command:
```
/process/(in)activate processName
```
     allows to activate/inactivate the processes one by one.

## ACTION INITALIZATION

   A newly introduced class, B1::ActionInitialization, instantiates and registers
   to Geant4 kernel all user action classes.

   While in sequential mode the action classes are instatiated just once,
   via invoking the method:
      B1::ActionInitialization::Build()
   in multi-threading mode the same method is invoked for each thread worker
   and so all user action classes are defined thread-local.

   A run action class is instantiated both thread-local
   and global that's why its instance is created also in the method
      B1::ActionInitialization::BuildForMaster()
   which is invoked only in multi-threading mode.

## PRIMARY GENERATOR

   The primary generator is defined in the B1::PrimaryGeneratorAction class.
   The default kinematics is a 6 MeV gamma, randomly distributed in front
   of the envelope across 80% of the transverse (X,Y) envelope size.
   This default setting can be changed via the Geant4 built-in commands
   of the G4ParticleGun class.

## DETECTOR RESPONSE

   This example demonstrates a simple scoring implemented directly
   in the user action classes.  Alternative ways of scoring via Geant4 classes
   can be found in the other examples.

   The energy deposited is collected step by step for a selected volume
   in B1::SteppingAction and accumulated event by event in B1::EventAction.

   At end of event, the value acummulated in B1::EventAction is added in B1::RunAction
   and summed over the whole run (see B1::EventAction::EndOfevent()).

   Total dose deposited is computed at B1::RunAction::EndOfRunAction(),
   and printed together with informations about the primary particle.
   In multi-threading mode the energy accumulated in G4Accumulable objects per
   workers is merged to the master in B1::RunAction::EndOfRunAction() and the final
   result is printed on the screen.

   G4Accumulable<G4double> type instead of G4double type is used for the B1::RunAction
   data members in order to facilitate merging of the values accumulated on workers
   to the master.  Currently the accumulables have to be registered to G4AccumulablesManager
   and G4AccumulablesManager::Merge() has to be called from the users code. This is planned
   to be further simplified with a closer integration of G4Accumulable classes in
   the Geant4 kernel next year.

   An example of creating and computing new units (e.g., dose) is also shown
   in the class constructor.

<hr>

## VISUALISATION

   The visualization manager is set via the G4VisExecutive class
   in the main() function in exampleB1.cc.
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

   The vis.mac macro in example B1 has additional commands
   that demonstrate additional functionality of the vis system, such as
   displaying text, axes, scales, date, logo and shows how to change
   viewpoint and style.  Consider copying these to other examples or
   your application.  To see even more commands use help or
   ls or browse the available UI commands in the Application
   Developers Guide, "Controlling Visualization from Commands".

## USER INTERFACES

   The user command interface is set via the G4UIExecutive class
   in the main() function in exampleB1.cc

   The selection of the user command interface is then done automatically
   according to the Geant4 configuration or it can be done explicitly via
   the third argument of the G4UIExecutive constructor (see exampleB4a.cc).

   The gui.mac macros are provided in examples B2, B4 and B5. This macro
   is automatically executed if Geant4 is built with any GUI session.
   It is also possible to customise the icons menu bar which is
   demonstrated in the icons.mac macro in example B5.

## HOW TO RUN

   - Execute exampleB1 in the 'interactive mode' with visualization
```
% exampleB1
and type in the commands from run1.mac line by line:
Idle> /control/verbose 2
Idle> /tracking/verbose 1
Idle> /run/beamOn 10
Idle> ...
Idle> exit
```
      or
```
Idle> /control/execute run1.mac
....
Idle> exit
```

   - Execute exampleB1 in the 'batch' mode from macro files
   (without visualization)
```
% exampleB1 run2.mac
% exampleB1 exampleB1.in > exampleB1.out
```


