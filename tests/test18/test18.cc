// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test18.cc,v 1.4 2000-06-14 12:33:01 flei Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - GGEVis.cc
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIGAG.hh"
#include "G4UIterminal.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif 

#include "Tst18GeometryConstruction.hh"
#include "Tst18PhysicsList.hh"
#include "Tst18EventAction.hh"
#include "Tst18RunAction.hh"
#include "Tst18SteppingAction.hh"
#include "Tst18PrimaryGeneratorAction.hh"
#include "Randomize.hh"

#include "g4std/vector"

G4bool drawEvent;
G4std::vector<G4String> Particles;
G4std::vector<G4double> Energies;
G4std::vector<G4double> Weights;
G4std::vector<G4double> Times;

#ifdef G4VIS_USE
#include "Tst18VisManager.hh"
#endif

int main(int argc,char** argv)
{

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes

  Tst18GeometryConstruction* Geometry = new Tst18GeometryConstruction;
  runManager->SetUserInitialization(Geometry);
  runManager->SetUserInitialization(new Tst18PhysicsList);

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new Tst18VisManager;
  visManager->Initialize();
#endif

  // set mandatory user action class
  runManager->SetUserAction(new Tst18PrimaryGeneratorAction);
  runManager->SetUserAction(new Tst18RunAction);
  runManager->SetUserAction(new Tst18EventAction);
  runManager->SetUserAction(new Tst18SteppingAction);
  // Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the User Interface manager 
    
      G4UImanager* UI = G4UImanager::GetUIpointer();  
      UI->ApplyCommand("/run/verbose 0");
      UI->ApplyCommand("/event/verbose 0");
      UI->ApplyCommand("/tracking/verbose 0");           
      UI->ApplyCommand("/process/verbose 0");           
      if (argc==1) {
	//     	G4UIsession * session = new G4UIGAG;
     	G4UIsession * session = new G4UIterminal;
      	session->SessionStart();
      	delete session;
  }
      else {
 
  //
  //
  // Create a pointer to the user interface manager.
  //
	//      G4UImanager *UI = G4UImanager::GetUIpointer();
      G4String command = "/control/execute ";
      for (int i=2; i<=argc; i++) {
      G4String macroFileName = argv[i-1];
      UI->ApplyCommand(command+macroFileName);
       }

  }                                  

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;

  return 0;
}








