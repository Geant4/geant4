///////////////////////////////////////////////////////////////////////////////
// File: HcalTestBeam96.cc
// Author: veronique.lefebure@cern.ch  11.1998
// 10/01 P.Arce  Reorganization of user actions   
// Description: Main function for Geant4 application HCAL Test-BEAM H2-96
///////////////////////////////////////////////////////////////////////////////

//#define VIS_USE

#include "HcalTestBeam96DetectorConstruction.hh"
#include "HcalTB96EndOfEventAction.hh"
#include "HcalTB96RunAction.hh"

#include "CMSPrimaryGeneratorAction.hh"
#include "HcalTestBeam99PhysicsList.hh"

#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#ifdef VIS_USE
  #include "CMSVisManager.hh"
#endif


int main(int argc,char** argv) {

#ifdef VIS_USE
  CMSVisManager *visManager = new CMSVisManager;
  visManager->Initialize();
#endif        

  G4RunManager * runManager = new G4RunManager;
  runManager->SetUserInitialization(new HcalTestBeam96DetectorConstruction);
  runManager->SetUserInitialization(new HcalTestBeam99PhysicsList);     

  ////////////////////////////
  //  User action classes.  //
  //  --------------------  //
  ////////////////////////////

  //////////////////////////////////
  // PRIMARY PARTICLEs GENERATION //
  //////////////////////////////////

  CMSPrimaryGeneratorAction* primaryGenerator = new CMSPrimaryGeneratorAction;
  runManager->SetUserAction(primaryGenerator);
  
  /////////
  // RUN //
  /////////

  runManager->SetUserAction(new HcalTB96RunAction);
  
  ///////////
  // EVENT //
  ///////////

  runManager->SetUserAction(new HcalTB96EndOfEventAction(primaryGenerator));
  
  G4UImanager * UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/OSCAR/generator/verbose 2");
  UI->ApplyCommand("/gun/position -1380. 0. 0. mm");
  UI->ApplyCommand("/gun/direction 1. 0. 0.");
  UI->ApplyCommand("/gun/energy 10 GeV");

  // Define (G)UI terminal for interactive mode
  if (argc==1) {
    G4UIsession * session = new G4UIterminal;
    cout <<" Run initializing ..."<<endl;
    UI->ApplyCommand("/process/verbose 0");
    UI->ApplyCommand("/run/verbose 2");
    UI->ApplyCommand("/run/initialize");
#ifdef VIS_USE
    //Ask for visualization, and some personal code.
    G4cout << "Which system do you want (0 for no visualization)? " << flush;
    G4String answer;
    cin >> answer;
    if (answer!="0") {
      G4String visSystem = "/vis/open ";
      visSystem+=answer;
      G4cout << endl << "Command: " << visSystem << endl;
      UI->ApplyCommand(visSystem);
    }
#endif
    UI->ApplyCommand("/tracking/storeTrajectory 0");
    cout <<"Now, please, apply beamOn command..."<<endl;
    session->SessionStart();    
    delete session;
  }
  else {
  // Batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

  delete runManager;

#ifdef VIS_USE
  delete visManager;
#endif

  return 0;
}







