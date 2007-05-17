//#define G4ANALYSIS_USE_ROOT 1

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyPhantomSD.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "HadrontherapyRunAction.hh"
#include "HadrontherapyMatrix.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "globals.hh"
#include "HadrontherapySteppingAction.hh"

#ifdef  G4ANALYSIS_USE
#include "HadrontherapyAnalysisManager.hh"
#endif

//#ifdef G4ANALYSIS_USE_ROOT
#include "rootAnalysis.hh"
//#endif

int main(int argc ,char ** argv)
{
  G4RunManager* runManager = new G4RunManager;

  // Initialize the geometry
  runManager -> SetUserInitialization(new HadrontherapyDetectorConstruction());
  
  // Initialize the physics 
  runManager -> SetUserInitialization(new HadrontherapyPhysicsList());
  
  // Initialize the primary particles  
  runManager -> SetUserAction(new HadrontherapyPrimaryGeneratorAction());

  // Initialize matrix 
  HadrontherapyMatrix* matrix = new HadrontherapyMatrix();
  matrix -> Initialize();


  // Optional UserActions: run, event, stepping
  HadrontherapyRunAction* pRunAction = new HadrontherapyRunAction();
  runManager -> SetUserAction(pRunAction);

  HadrontherapyEventAction* pEventAction = new HadrontherapyEventAction(matrix);
  runManager -> SetUserAction(pEventAction);


  HadrontherapySteppingAction* pSteppingAction = new HadrontherapySteppingAction(pRunAction); 
  runManager -> SetUserAction(pSteppingAction);    


#ifdef G4ANALYSIS_USE
  //  HadrontherapyAnalysisManager* analysis = 
  //  HadrontherapyAnalysisManager::getInstance();
  //  analysis->book();
#endif

  //#ifdef G4ANALYSIS_USE_ROOT
  //  rootAnalysis* analyseWithRoot = new rootAnalysis;
    rootAnalysis* analysis = rootAnalysis::getInstance();
  //  rootAnalysis* analyseWithRoot = rootAnalysis::getInstance();
  //  analyseWithRoot->book();
  analysis->book();
  #//endif
  
#ifdef G4VIS_USE
  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  G4UIsession* session = 0;
  if (argc == 1)   // Define UI session for interactive mode.
    {
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif

      //      session = new G4UIterminal();
    } 

  // Get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  if (session)   // Define UI session for interactive mode.
    { 
      G4cout<<" UI session starts ..."<< G4endl;
      //      UI -> ApplyCommand("/control/execute defaultMacro.mac");    
      UI -> ApplyCommand("/control/execute hadronics.mac");    
      session -> SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI -> ApplyCommand(command + fileName);
    }  

  matrix -> TotalEnergyDeposit();

#ifdef G4ANALYSIS_USE
  //  analysis->finish();
#endif


  //#ifdef G4ANALYSIS_USE_ROOT
  //   analyseWithRoot->finish();
       analysis->finish();
   //#endif
  
  // Job termination
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;

  return 0;
}
