#include "B01VisApplication.hh"
#include "B01VSimulation.hh"
#include "G4VisManager.hh"
#include "B01VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

B01VisApplication::B01VisApplication()
{}
B01VisApplication::~B01VisApplication(){}

void B01VisApplication::RunSimulation(B01VSimulation *sim){
  G4UIsession *session = CreateSession();


  sim->PrepareSampling();

  fVisRun.SetDetector(sim->GetMassGeometry());
  fVisRun.Initialize();
  
  sim->ConfigureSampling();

  G4UImanager::GetUIpointer()->ApplyCommand("/control/execute vis.mac");
  session->SessionStart();
  
  sim->PostRun(&G4std::G4cout);
}





G4UIsession *B01VisApplication::CreateSession(){
  G4std::G4cout << "========== createing B01VisManager" << G4endl;
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new B01VisManager;
  visManager->Initialize();



  G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
  session = new G4UIterminal(new G4UItcsh);      
#else
  session = new G4UIterminal();
#endif    

  return session;
}



