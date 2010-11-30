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
//
// $Id: cell.cc,v 1.32 2006/09/20 15:30:15 guatelli
//
// 

#include "globals.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "CellDetectorConstruction.hh"
#include "CellPhysicsList.hh"
#include "CellPrimaryGeneratorAction.hh"
#include "CellRunAction.hh"
#include "CellEventAction.hh"
#include "CellSteppingAction.hh"
#include "CellSteppingVerbose.hh"
#include "CellAnalysisManager.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv) {
 	
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);   
  G4int seed = time(0);
  CLHEP::HepRandom ::setTheSeed(seed);
 
  G4RunManager * runManager = new G4RunManager;

  CellDetectorConstruction* cellDetector = new CellDetectorConstruction();
  runManager->SetUserInitialization(cellDetector);  

  CellPhysicsList* cellPhysics = new CellPhysicsList();
  runManager->SetUserInitialization(cellPhysics);

#ifdef G4VIS_USE 
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  CellPrimaryGeneratorAction* cellPrimaryParticle = new CellPrimaryGeneratorAction(); 
  runManager->SetUserAction(cellPrimaryParticle);

  CellRunAction* cellRun=new CellRunAction(); 
  runManager->SetUserAction(cellRun);  

  CellEventAction *cellEventAction = new CellEventAction();
 
  runManager->SetUserAction(cellEventAction);
     
  CellSteppingAction* cellSteppingAction = new CellSteppingAction(cellPrimaryParticle,
                                                                cellRun,
                                                                cellDetector);
  runManager->SetUserAction(cellSteppingAction);

  CellAnalysisManager* analysis = CellAnalysisManager::getInstance();
 
  // Set the name of the generated .hbk file
  analysis -> book("cell");
  
  //get the pointer to the User Interface manager 
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  
  
  if (argc == 1)
    // Define (G)UI terminal for interactive mode  
    { 
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      UImanager->ApplyCommand("/control/execute default.mac");
      ui->SessionStart();
      delete ui;
#endif
    }
  else
    // Batch mode
    {     
      G4String command =("/control/execute ");
      G4String fileName = argv[1];
      G4cout <<"macro --> "<< fileName << G4endl;
      UImanager->ApplyCommand(command+fileName);
    }

  analysis->finish();

#ifdef G4VIS_USE
  delete visManager;
#endif

   delete runManager;

  return 0;
}
