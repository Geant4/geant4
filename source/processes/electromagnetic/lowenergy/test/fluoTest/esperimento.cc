#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif
#include "Randomize.hh"
#ifdef G4VIS_USE
#include "myVisManager.hh"
#endif
#include "myDetectorConstruction.hh"
#include "myPhysicsList.hh"
#include "myPrimaryGeneratorAction.hh"
#include "myRunAction.hh"
#include "myEventAction.hh"
#include "mySteppingAction.hh"
#include "mySteppingVerbose.hh"

#ifdef G4ANALYSIS_USE
#include "myAnalysisManager.hh"
#endif

/* This global file is used to store relevant data for
   analysis with external tools */
G4std::ofstream outFile;


int main(int argc,char** argv) {

  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);

  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new mySteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  myDetectorConstruction* detector = new myDetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new myPhysicsList);
  
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
  //visualization manager
   G4VisManager* visManager = new myVisManager;
  visManager->Initialize();
#endif
  
#ifdef G4ANALYSIS_USE
  // Creation of the analysis manager
  myAnalysisManager* analysisMgr = new myAnalysisManager(detector);
#endif
  
 // Set optional user action classes
#ifdef G4ANALYSIS_USE
  myEventAction* eventAction = 
    new myEventAction(analysisMgr);
  myRunAction* runAction =
    new myRunAction(analysisMgr);
  mySteppingAction* stepAction = 
    new mySteppingAction(detector,analysisMgr); 
 //#else 
  //myEventAction* eventAction = new myEventAction();
  //myRunAction* runAction = new myRunAction();
#endif 

// set user action classes
  runManager->SetUserAction(new myPrimaryGeneratorAction(detector));
  
  runManager->SetUserAction(eventAction); 
  runManager->SetUserAction(runAction);
#ifdef G4ANALYSIS_USE 
runManager->SetUserAction(stepAction);
#else
runManager->SetUserAction(new mySteppingAction);  
#endif
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
#ifdef G4ANALYSIS_USE
  delete analysisMgr;  
 #endif
delete runManager;
  return 0;
}

 

