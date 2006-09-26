//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include <stdexcept>

#include "globals.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4TransportationManager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "G4HumanPhantomConstruction.hh"
#include "G4HumanPhantomPhysicsList.hh"
#include "G4HumanPhantomPrimaryGeneratorAction.hh"
#include "G4HumanPhantomSteppingAction.hh"
#include "G4HumanPhantomEventAction.hh"
#include "G4HumanPhantomRunAction.hh"
#include "G4HumanPhantomEnergyDeposit.hh"
#include  "G4HumanPhantomAnalysisManager.hh"

int main(int argc,char** argv)
{
   G4cout << " starting run " << G4endl;

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;
  
  // set mandatory initialization classes
  G4HumanPhantomConstruction* userPhantom = new G4HumanPhantomConstruction;
  runManager->SetUserInitialization(userPhantom);
  G4cout << " useraction geometry " << G4endl;


  runManager->SetUserInitialization(new G4HumanPhantomPhysicsList);
  G4cout << " useraction PhysicsList " << G4endl;

  runManager->SetUserAction(new G4HumanPhantomPrimaryGeneratorAction);
   G4cout << " useraction Primary Generator " << G4endl;

  G4UIsession* session=0;

  if (argc==1)   
  // Define UI session for interactive mode. 
    { 
     // G4UIterminal is a (dumb) terminal
#ifdef G4UI_USE_TCSH
     session = new G4UIterminal(new G4UItcsh);
#else
     session = new G4UIterminal();
#endif
     }

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  runManager->Initialize();
  G4cout << "runManager->Initialize();" << G4endl;

  // set mandatory user action class
  
  runManager->SetUserAction(new G4HumanPhantomRunAction);

  G4HumanPhantomEnergyDeposit* energyTotal = new G4HumanPhantomEnergyDeposit();

  G4HumanPhantomEventAction *eventAction = new G4HumanPhantomEventAction(energyTotal);
  runManager->SetUserAction(eventAction);

  runManager->SetUserAction(new G4HumanPhantomSteppingAction(eventAction)); 


  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (session)   // Define UI session for interactive mode.
    { 
      G4cout << " UI session starts ..." << G4endl;
      UI -> ApplyCommand("/control/execute run.mac");    
      session -> SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI -> ApplyCommand(command+fileName);
    }     

   energyTotal->TotalEnergyDeposit();

#ifdef G4ANALYSIS_USE
  G4HumanPhantomAnalysisManager* analysis = G4HumanPhantomAnalysisManager::getInstance();
  if (analysis)analysis->finish();
#endif

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;
  return 0;
}
