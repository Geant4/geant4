// $Id: README.TXT,v 1.1 1999-01-07 16:14:15 gunter Exp $

 GEANT4
 
 This directory contains several examples of GEANT4 main programs,
 from very simple cases (a purely batch program), 
 to interactive programs based on the G4 command line interface,
 to more complex cases programs touching most of the G4 functionality.
 
 Those main programs provide just examples of initialization and usage
 of the GEANT4 toolkit classes, but the user is also free to define his
 own way and functions to initialize G4 and write the main program.

 A typical GEANT4 main program looks like the following:

int main() {

  // Create Run Manager  
  G4RunManager * runManager = new G4RunManager;

  // Register User Classes to the RunManager
    // Mandatory classes  -----------------------
      // Detector geometry
      runManager-> set_userInitialization(new MyDetectorConstruction);
      // Physics List
      runManager-> set_userInitialization(new MyPhysicsList);
      // Primary Generator
      runManager->set_userAction(new MyPrimaryGeneratorAction);
  
    // Optiolnal classes  -----------------------
      // User Actions
      runManager->set_userAction(new MyRunAction);
      runManager->set_userAction(new MyEventAction);
      runManager->set_userAction(new MyStackingAction);
      runManager->set_userAction(new MyTrackingAction);
      runManager->set_userAction(new MySteppingAction);


  // Define (G)UI terminal for interactive mode
  G4UIsession * session = new G4UIterminal;

  // User interactions
  session->sessionStart();
  
  // Termination
  delete session;
  delete runManager;
  return 0;
}
 
 First of all, user must create the RunManager. The RunManager controles 
 run sequence by receiving messages from the user via UIsession.
 The following are a part of commands for run controle; 
    /run/initialize *          Initialize G4 kernel.
    /run/beamOn *              Start a Run.
    /run/verbose *             Set the verbose level of G4RunManager.
    /run/abort *               Abort current run processing.
 For example, event loop will start by using "/run/beamOn" commands.

 In order to execute simulation, user must provide geometrical 
 configuration of his own detector. In addition information of 
 primary events must be given with a list of particle types and
 processes for them.
 
 So, user must provide his own classes derived from the following
 three abstract classes by implementing their pure virtual functions
 and register those user classes to the RunManager. 

 G4VUserDetectorConstruction       - Detector Geometry, Materials
   pure virtual functions
     G4VPhysicalVolume* construct() 
	   - construct detectors;
             normally implemented to serve as the entry point 
             for the tree of methods describing solids, volumes, 
             materials and sensitive detectors.

 G4VUserPhysicsList                - Particle types and Processes 
   pure virtual functions
     void constructParticle() 
           - construct particles;
             normally implemented to select desired particle types.  
     void constructPhysics() 
           - construct procesess;
	     normally implemented to select desired physics processes
             for each particle types and register them to the ProcessManager.
     void setCuts(G4double aValue)
	   - sets a cut value; 
	     normally implemented to set cut value in range for all 
             particle types.
 
 G4VUserPrimaryGeneratorAction     - Event Generator selection
   pure virtual functions
     void generatePrimaries(G4Event* anEvent)
           - generate a event with primary particles
	     normally implemented to select the desired event generation 
             mechanism, such as the ParticleGun or the PYTHIA interface.

 In addition to the above mandatory classes, there are 5 user classes
 to customize the default functionality of GEANT4 simulation. 

 G4UserRunAction                   - Actions for each Run
 G4UserEventAction                 - Actions for each Event
 G4UserStackingAction              - Tracks Stacking selection
 G4UserTrackingAction              - Actions for each Track
 G4UserSteppingAction              - Actions for each Step
 
 The virtual functions belonging to the classes above are:

 G4UserRunAction                   - void beginOfRunAction(G4Run*)
                                   - void endOfRunAction(G4Run*)
 
 G4UserEventAction                 - void beginOfEventAction(G4Event*)
                                   - void endOfEventAction(G4Event*)
 
 G4UserStackingAction              - G4ClassificationOfNewTrack
                                         classifyNewTrack(G4Track *const)
                                   - void newStage()
                                   - void prepareNewEvent()
 
 G4UserTrackingAction              - void preUserTrackingAction()
                                   - void postUserTrackingAction()
  
 G4UserSteppingAction              - void userSteppingAction()
 
 Finally, more details can be found in the header files and in the source
 code relative to the classes outlined above.






