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

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
#include "HadrontherapyVisManager.hh"
#endif

#include "HadrontherapyEventAction.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyPhantomSD.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "G4SDManager.hh"
#include "HadrontherapyRunAction.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"

G4String sensitiveDetectorName;
G4int k;
G4int i;
G4int j;
G4double matrix[40][40][40]; // dimansions of the output matrix

// ----------------------------------------------------------------
int main(int argc ,char ** argv)

{
HepRandom::setTheEngine(new RanecuEngine);
G4int seed = time(0);
HepRandom :: setTheSeed(seed);

G4RunManager* pRunManager = new G4RunManager;

G4String sensitiveDetectorName = "Phantom";

HadrontherapyDetectorConstruction  *pDetectorConstruction = 
  new  HadrontherapyDetectorConstruction(sensitiveDetectorName);


pRunManager->SetUserInitialization(pDetectorConstruction);
pRunManager->SetUserInitialization(new HadrontherapyPhysicsList(pDetectorConstruction));

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new HadrontherapyVisManager;
  visManager->Initialize();
#endif
  
  
G4UIsession* session = 0;
 if (argc == 1)   // Define UI session for interactive mode.
   {
     session = new G4UIterminal();
   }

  
pRunManager->SetUserAction(new HadrontherapyRunAction(sensitiveDetectorName));

HadrontherapyEventAction *pEventAction = new HadrontherapyEventAction();
pRunManager->SetUserAction(pEventAction );

HadrontherapyPrimaryGeneratorAction* primary = 
  new HadrontherapyPrimaryGeneratorAction(pDetectorConstruction);
pRunManager->SetUserAction(primary);


   //Initialize G4 kernel
pRunManager->Initialize();

  // get the pointer to the User Interface manager 
G4UImanager* UI = G4UImanager::GetUIpointer();  
 if (session)   // Define UI session for interactive mode.
   { 
     G4cout<<" UI session starts ..."<< G4endl;
     UI->ApplyCommand("/control/execute VisualisationMacro.mac");    
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

  delete pRunManager;

  return 0;
}
