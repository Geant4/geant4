//  XrayTel.cc

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIGAG.hh"
#include "G4UIXm.hh" 

#include "XrayTelDetectorConstruction.hh"
#include "XrayTelPhysicsList.hh"
#include "XrayTelVisManager.hh"
#include "XrayTelEventAction.hh"
#include "XrayTelRunAction.hh"
#include "XrayTelSteppingAction.hh"
#include "XrayTelPrimaryGeneratorAction.hh"

#include <iostream.h>

#include "g4std/vector"

G4bool drawEvent;
G4std::vector<G4String> EnteringParticles;
G4std::vector<G4double> EnteringEnergy;
G4std::vector<G4ThreeVector> EnteringDirection;

int main( int argc, char** argv )
{
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  XrayTelDetectorConstruction* telescope = new XrayTelDetectorConstruction;
  runManager->SetUserInitialization( telescope );
  runManager->SetUserInitialization(new XrayTelPhysicsList);

  // set mandatory user action class
  runManager->SetUserAction(new XrayTelPrimaryGeneratorAction( telescope ));
  runManager->SetUserAction(new XrayTelRunAction);
  runManager->SetUserAction(new XrayTelEventAction);
  runManager->SetUserAction(new XrayTelSteppingAction);
   
  // visualization manager
  G4VisManager* visManager = new XrayTelVisManager;
  visManager->Initialize();    

  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager *UI = G4UImanager::GetUIpointer();  
  if ( argc==1 ){
    G4UIsession * session = new G4UIGAG;
    session->SessionStart();
    delete session;
  }
  else {
    // Create a pointer to the User Interface manager 
    G4String command = "/control/execute ";
    for (int i=2; i<=argc; i++) {
       G4String macroFileName = argv[i-1];
       UI->ApplyCommand(command+macroFileName);
    }
  }                                  

  // job termination
  delete visManager;
  delete runManager;
  return 0;
}

