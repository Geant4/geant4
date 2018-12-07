-------------------------------------------------------------------

     =========================================================
     Geant4 - an Object-Oriented Toolkit for Simulation in HEP
     =========================================================

                            Example B1
                            -----------

 This example demonstrates a very simple application where an energy
 deposit is accounted in user actions and their associated objects
 and a dose in a selected volume is calculated. 

	
 1- GEOMETRY DEFINITION
	
   The geometry is constructed in the B1DetectorConstruction class.
   The setup consists of a an envelope of box shape containing two
   volumes: a spherical cone and a trapezoid.

   In this example we use  some common materials materials for medical
   applications. The envelope is made of water and the two inner volumes
   are made from tissue and bone materials.  
   The materials are created with the help of the G4NistManager class,
   which allows to build a material from the NIST database using their
   names. All available materials can be found in the Geant4 User's Guide
   for Application Developers, Appendix 10: Geant4 Materials Database.
		
 2- PHYSICS LIST
 
   The particle's type and the physic processes which will be available
   in this example are set in the QBBC physics list. This physics list 
   requires data files for electromagnetic and hadronic processes.
   See more on installation of the datasets in Geant4 Installation Guide,
   Chapter 3.3: Note On Geant4 Datasets:
   http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides
                                           /InstallationGuide/html/ch03s03.html
   The following datasets: G4LEDATA, G4LEVELGAMMADATA, G4NEUTRONXSDATA,
   G4SAIDXSDATA and G4ENSDFSTATEDATA are mandatory for this example.

   In addition the build-in interactive command:
               /process/(in)activate processName
   allows to activate/inactivate the processes one by one.
   
 3- ACTION INITALIZATION

   A newly introduced class, B1ActionInitialization, instantiates and registers 
   to Geant4 kernel all user action classes.

   While in sequential mode the action classes are instatiated just once,
   via invoking the method:
      B1ActionInitialization::Build() 
   in multi-threading mode the same method is invoked for each thread worker
   and so all user action classes are defined thread-local.

   A run action class is instantiated both thread-local 
   and global that's why its instance is created also in the method
      B1ActionInitialization::BuildForMaster() 
   which is invoked only in multi-threading mode.
  	 
 4- PRIMARY GENERATOR
  
   The primary generator is defined in the B1PrimaryGeneratorAction class.
   The default kinematics is a 6 MeV gamma, randomly distributed in front
   of the envelope across 80% of the transverse (X,Y) envelope size. 
   This default setting can be changed via the Geant4 built-in commands 
   of the G4ParticleGun class.
     
 5- DETECTOR RESPONSE

   This example demonstrates a simple scoring implemented directly
   in the user action classes.  Alternative ways of scoring via Geant4 classes 
   can be found in the other examples.
   
   The energy deposited is collected step by step for a selected volume
   in B1SteppingAction and accumulated event by event in B1EventAction.

   At end of event, the value acummulated in B1EventAction is added in B1RunAction
   and summed over the whole run (see B1EventAction::EndOfevent()).

   Total dose deposited is computed at B1RunAction::EndOfRunAction(), 
   and printed together with informations about the primary particle.
   In multi-threading mode the energy accumulated in G4Accumulable objects per
   workers is merged to the master in B1RunAction::EndOfRunAction() and the final
   result is printed on the screen.
   
   G4Parameter<G4double> type instead of G4double type is used for the B1RunAction
   data members in order to facilitate merging of the values accumulated on workers 
   to the master.  Currently the accumulables have to be registered to G4AccumulablesManager
   and G4AccumulablesManager::Merge() has to be called from the users code. This is planned
   to be further simplified with a closer integration of G4Accumulable classes in
   the Geant4 kernel next year.

   An example of creating and computing new units (e.g., dose) is also shown 
   in the class constructor. 

 The following paragraphs are common to all basic examples

 A- VISUALISATION

   The visualization manager is set via the G4VisExecutive class
   in the main() function in exampleB1.cc.    
   The initialisation of the drawing is done via a set of /vis/ commands
   in the macro vis.mac. This macro is automatically read from
   the main function when the example is used in interactive running mode.

   By default, vis.mac opens an OpenGL viewer (/vis/open OGL).
   The user can change the initial viewer by commenting out this line
   and instead uncommenting one of the other /vis/open statements, such as
   HepRepFile or DAWNFILE (which produce files that can be viewed with the
   HepRApp and DAWN viewers, respectively).  Note that one can always
   open new viewers at any time from the command line.  For example, if
   you already have a view in, say, an OpenGL window with a name
   "viewer-0", then
      /vis/open DAWNFILE
   then to get the same view
      /vis/viewer/copyView viewer-0
   or to get the same view *plus* scene-modifications
      /vis/viewer/set/all viewer-0
   then to see the result
      /vis/viewer/flush

   The DAWNFILE, HepRepFile drivers are always available
   (since they require no external libraries), but the OGL driver requires
   that the Geant4 libraries have been built with the OpenGL option.

   From Release 9.6 the vis.mac macro in example B1 has additional commands
   that demonstrate additional functionality of the vis system, such as
   displaying text, axes, scales, date, logo and shows how to change
   viewpoint and style.  Consider copying these to other examples or
   your application.  To see even more commands use help or
   ls or browse the available UI commands in the Application
   Developers Guide, Section 7.1.

   For more information on visualization, including information on how to
   install and run DAWN, OpenGL and HepRApp, see the visualization tutorials,
   for example,
   http://geant4.slac.stanford.edu/Presentations/vis/G4[VIS]Tutorial/G4[VIS]Tutorial.html
   (where [VIS] can be replaced by DAWN, OpenGL and HepRApp)

   The tracks are automatically drawn at the end of each event, accumulated
   for all events and erased at the beginning of the next run.

 B- USER INTERFACES
 
   The user command interface is set via the G4UIExecutive class
   in the main() function in exampleB1.cc 
   The selection of the user command interface is then done automatically 
   according to the Geant4 configuration or it can be done explicitly via 
   the third argument of the G4UIExecutive constructor (see exampleB4a.cc). 
 
 C- HOW TO RUN

    - Execute exampleB1 in the 'interactive mode' with visualization:
        % ./exampleB1
      and type in the commands from run1.mac line by line:  
        Idle> /control/verbose 2
        Idle> /tracking/verbose 1
        Idle> /run/beamOn 10 
        Idle> ...
        Idle> exit
      or
        Idle> /control/execute run1.mac
        ....
        Idle> exit

    - Execute exampleB1  in the 'batch' mode from macro files 
      (without visualization)
        % ./exampleB1 run2.mac
        % ./exampleB1 exampleB1.in > exampleB1.out

	
