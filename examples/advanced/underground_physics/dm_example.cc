// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// 
// --------------------------------------------------------------
//   GEANT 4 - Project 2001 - Underground Dark Matter Detector
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "ZIIIVisManager.hh"
#endif

#include "ZIIIDetectorConstruction.hh"
#include "ZIIIPhysicsList.hh"
#include "ZIIIPrimaryGeneratorAction.hh"
#include "ZIIIRunAction.hh"
#include "ZIIIEventAction.hh"
#include "ZIIISteppingAction.hh"
#include "ZIIISteppingVerbose.hh"

#include "g4std/vector"
//G4String filename;
G4bool drawEvent;
G4std::vector<G4String> Particles;
G4std::vector<G4double> Energies;
G4std::vector<G4double> Weights;
G4std::vector<G4double> Times;


int main(int argc,char** argv) {

  //  G4cout << " The results file name = "<<G4endl ;
  //G4cin >> filename;


  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new ZIIISteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  ZIIIDetectorConstruction* detector = new ZIIIDetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new ZIIIPhysicsList);
  
 G4UIsession* session=0;
  
  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_XM
      session = new G4UIXm(argc,argv);
#else           
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif
#endif
    }
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new ZIIIVisManager;
  visManager->Initialize();
#endif
    
  // set user action classes
  runManager->SetUserAction(new ZIIIPrimaryGeneratorAction(detector));
  runManager->SetUserAction(new ZIIIRunAction);
  runManager->SetUserAction(new ZIIIEventAction);
  runManager->SetUserAction(new ZIIISteppingAction);
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (session)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      UI->ApplyCommand("/control/execute initInter.mac");    
#ifdef G4UI_USE_XM
      // Customize the G4UIXm menubar with a macro file :
      UI->ApplyCommand("/control/execute gui.mac");
#endif
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

