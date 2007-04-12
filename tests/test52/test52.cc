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
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// $Id: test52.cc,v 1.1 2007-04-12 12:00:40 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "globals.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Tst52DetectorConstruction.hh"
#include "Tst52PhysicsList.hh"
#include "Tst52PrimaryGeneratorAction.hh"
#include "Tst52RunAction.hh"
#include "Tst52EventAction.hh"
#include "Tst52SteppingAction.hh"
#include "Tst52SteppingVerbose.hh"
#include "Tst52AnalysisManager.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

int main(int argc,char** argv) {
 	
  // CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);   
  // G4int seed = time(0);
  // CLHEP::HepRandom ::setTheSeed(seed);
 
  G4RunManager * runManager = new G4RunManager;

  Tst52DetectorConstruction* tst52Detector = new Tst52DetectorConstruction();
  runManager->SetUserInitialization(tst52Detector);  

  Tst52PhysicsList* tst52Physics = new Tst52PhysicsList();
  runManager->SetUserInitialization(tst52Physics);

#ifdef G4VIS_USE 
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  Tst52PrimaryGeneratorAction* tst52PrimaryParticle = new Tst52PrimaryGeneratorAction(); 
  runManager->SetUserAction(tst52PrimaryParticle);

  Tst52RunAction* tst52Run=new Tst52RunAction(); 
  runManager->SetUserAction(tst52Run);  

  Tst52EventAction *tst52EventAction = new Tst52EventAction();
 
  runManager->SetUserAction(tst52EventAction);
     
  Tst52SteppingAction* tst52SteppingAction = new Tst52SteppingAction(tst52PrimaryParticle,
                                                                tst52Run,
                                                                tst52Detector);
  runManager->SetUserAction(tst52SteppingAction);

  Tst52AnalysisManager* analysis = Tst52AnalysisManager::getInstance();
 
  analysis->book();
  
  //get the pointer to the User Interface manager 
  G4UImanager * UI = G4UImanager::GetUIpointer();  
  
  if (argc == 1)
    // Define (G)UI terminal for interactive mode  
    { 
      // G4UIterminal is a (dumb) terminal.
      G4UIsession * session = 0;

#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    
      UI->ApplyCommand("/control/execute default.mac");
      session->SessionStart();
      delete session;
    }
  else
    // Batch mode
    {     
      G4String command =("/control/execute ");
      G4String fileName = argv[1];
      G4cout <<"macro --> "<< fileName << G4endl;
      UI->ApplyCommand(command+fileName);
    }

  analysis->finish();

#ifdef G4VIS_USE
  delete visManager;
#endif

   delete runManager;

  return 0;
}
