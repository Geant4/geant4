#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif
#include "Randomize.hh"
#ifdef G4VIS_USE
#include "FluoTestVisManager.hh"
#endif
#include "FluoTestDetectorConstruction.hh"
#include "FluoTestPhysicsList.hh"
#include "FluoTestPrimaryGeneratorAction.hh"
#include "FluoTestRunAction.hh"
#include "FluoTestEventAction.hh"
#include "FluoTestSteppingAction.hh"
#include "FluoTestSteppingVerbose.hh"

#ifdef G4ANALYSIS_USE
#include "FluoTestAnalysisManager.hh"
#endif

int main(int argc,char** argv) {

  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);

  //FluoTest Verbose output class
  G4VSteppingVerbose::SetInstance(new FluoTestSteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  FluoTestDetectorConstruction* detector = new FluoTestDetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new FluoTestPhysicsList);
  
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
   G4VisManager* visManager = new FluoTestVisManager;
  visManager->Initialize();
#endif
  
#ifdef G4ANALYSIS_USE
  // Creation of the analysis manager
  FluoTestAnalysisManager* analysisMgr = new FluoTestAnalysisManager(detector);
#endif
  
 // Set optional user action classes
#ifdef G4ANALYSIS_USE
  FluoTestEventAction* eventAction = 
    new FluoTestEventAction(analysisMgr);
  FluoTestRunAction* runAction =
    new FluoTestRunAction(analysisMgr);
  // FluoTestSteppingAction* stepAction = 
  // new FluoTestSteppingAction(detector,analysisMgr);
 FluoTestSteppingAction* stepAction = 
  new FluoTestSteppingAction(analysisMgr);
 #else 
 FluoTestEventAction* eventAction = new FluoTestEventAction();
  FluoTestRunAction* runAction = new FluoTestRunAction();
  FluoTestSteppingAction* stepAction = new FluoTestSteppingAction();
#endif 

// set user action classes
  runManager->SetUserAction(new FluoTestPrimaryGeneratorAction(detector));
  
  runManager->SetUserAction(eventAction); 
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(stepAction);

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
  /*
#ifdef G4VIS_USE
   delete visManager;
  #endif
  */
#ifdef G4ANALYSIS_USE
  delete analysisMgr;  
 #endif
delete runManager;
  return 0;
}

 


