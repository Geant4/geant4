//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Authors: S.Guatelli, Matthew Large and A. Malaroda, University of Wollongong
// 

#include "G4UImanager.hh"
#include "G4UIsession.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "ICRP110PhantomConstruction.hh"
#include "ICRP110PhantomActionInitialization.hh"
#include "G4ScoringManager.hh"
#include "ICRP110UserScoreWriter.hh"
#include "ICRP110PhantomVisAction.hh"
#include "QGSP_BIC_HP.hh"
#include "G4RunManagerFactory.hh"

int main(int argc,char** argv)
{
 auto* runManager = G4RunManagerFactory::CreateRunManager();
 G4int nThreads = 4;
 runManager->SetNumberOfThreads(nThreads);
  
  // Activate UI-command base scorer
 G4ScoringManager* scorerManager = G4ScoringManager::GetScoringManager();
 scorerManager -> SetVerboseLevel(1);

//====================================================================
// Un-comment this line for user defined score writer
    scorerManager -> SetScoreWriter(new ICRP110UserScoreWriter());
//====================================================================

 // Set mandatory initialization classes
  auto userPhantom = new ICRP110PhantomConstruction();
  runManager -> SetUserInitialization(userPhantom);
  
  runManager -> SetUserInitialization(new QGSP_BIC_HP());

 // runManager -> SetUserInitialization(new ICRP110PhantomPhysicsList);

  G4VisManager* visManager = new G4VisExecutive;
  visManager->RegisterRunDurationUserVisAction
  ("phantom",new ICRP110PhantomVisAction(userPhantom));
  visManager -> Initialize();
 
  auto actions = new ICRP110PhantomActionInitialization();
  runManager -> SetUserInitialization(actions);

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc==1)   // Define UI session for interactive mode.
    { 
      G4cout << " UI session starts ..." << G4endl;
      auto ui = new G4UIExecutive(argc, argv);
      UImanager -> ApplyCommand("/control/execute vis.mac");     
      ui -> SessionStart();
      delete ui;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager -> ApplyCommand(command+fileName);
    }     

delete visManager;

delete runManager;

return 0;
}
