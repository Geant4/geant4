#include "MLRunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
#include "MLVisManager.hh"
#endif

#include "MLGeometryConstruction.hh"
#include "MLPhysicsList.hh"
#include "MLPrimaryGeneratorAction.hh"
#include "MLAnalysisManager.hh"
#include "MLRunAction.hh"
#include "MLEventAction.hh"
#include "MLSteppingAction.hh"
#include "MLGeometryType.hh"

int main(int argc, char** argv)
{
  HepRandom::setTheEngine(new RanecuEngine);

  // Construct the default run manager
  MLRunManager* runManager = new MLRunManager;
  
  // Set mandatory user initialization classes
  MLGeometryConstruction* geometry = new MLGeometryConstruction;
  runManager->SetUserInitialization(geometry);
  runManager->SetUserInitialization(new MLPhysicsList);

  // Creation of the analysis manager
  MLAnalysisManager* analysis = MLAnalysisManager::getInstance();

    
  // Set optional user action classes
  MLEventAction* eventAction = new MLEventAction(geometry);
  MLRunAction* runAction = new MLRunAction();
  
  // Set mandatory user action classes
  runManager->SetUserAction(new MLPrimaryGeneratorAction());
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(new MLSteppingAction(geometry));
  
#ifdef G4VIS_USE
  // Visualization manager
  G4VisManager* visManager = new MLVisManager;
  visManager->Initialize();
#endif
  
  // Initialize G4 kernel
  // runManager->Initialize();
  
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  G4double xpos = -(geometry->GetWorldSizeZ())/2. - 5e-4;
  if ( geometry->GetShape() == SPHERE ) xpos *= 2 ;
  G4String command = "/gps/centre ";
  char x[10];
  std::ostrstream os(x,10);
  os << xpos << '\0'; 
  G4String xs = x;
  if ( geometry->GetShape() == SPHERE ) {
    UI->ApplyCommand(command+xs+" 0. 0. mm");
    UI->ApplyCommand("/gps/type Point");
    UI->ApplyCommand("/gps/angrot1 0. 0. 1.");
    UI->ApplyCommand("/gps/angrot2 0. 1. 0.");
    UI->ApplyCommand("/gps/particle proton");
    UI->ApplyCommand("/gps/direction 1 0 0");
    UI->ApplyCommand("/gps/energy 100 MeV");
  }else{
    UI->ApplyCommand(command+"0. 0. "+xs+" mm");
    UI->ApplyCommand("/gps/type Point");
    UI->ApplyCommand("/gps/angrot1 1. 0. 0.");
    UI->ApplyCommand("/gps/angrot2 0. -1. 0.");
    UI->ApplyCommand("/gps/particle proton");
    UI->ApplyCommand("/gps/direction 0 0 1");
    UI->ApplyCommand("/gps/energy 100 MeV");
  }

  if (argc==1)   // Define UI session for interactive mode.
    {
      G4UIsession* session=0;
#ifdef G4UI_USE_XM
      session = new G4UIXm(argc,argv);
#else           
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);
#else
      session = new G4UIterminal();
#endif
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

  // Job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete analysis;
  delete runManager;

  return 0;
}








