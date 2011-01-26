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
//
// $Id: test50.cc,v 1.34 2010-06-07 10:08:39 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "globals.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Tst50DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "Tst50PrimaryGeneratorAction.hh"
#include "Tst50RunAction.hh"
#include "Tst50EventAction.hh"
#include "Tst50SteppingAction.hh"
#include "Tst50SteppingVerbose.hh"
#include "Tst50AnalysisManager.hh"
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

  Tst50DetectorConstruction* tst50Detector = new Tst50DetectorConstruction();
  runManager->SetUserInitialization(tst50Detector);  

  PhysicsList* tst50Physics = new PhysicsList();
  runManager->SetUserInitialization(tst50Physics);

  Tst50PrimaryGeneratorAction* tst50PrimaryParticle = new Tst50PrimaryGeneratorAction(); 
  runManager->SetUserAction(tst50PrimaryParticle);

  Tst50RunAction* tst50Run=new Tst50RunAction(); 
  runManager->SetUserAction(tst50Run);  

  Tst50EventAction *tst50EventAction = new Tst50EventAction();
 
  runManager->SetUserAction(tst50EventAction);
     
  Tst50SteppingAction* tst50SteppingAction 
    = new Tst50SteppingAction(tst50PrimaryParticle, tst50Run, tst50Detector);
  runManager->SetUserAction(tst50SteppingAction);

#ifdef G4ANALYSIS_USE
  Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
 
  if (argc == 1) analysis->book("test50");
  else           analysis->book(argv[1]);
#endif
  
  //get the pointer to the User Interface manager 
  G4UImanager * UI = G4UImanager::GetUIpointer();  
  
  if (argc == 1)
    // Define (G)UI terminal for interactive mode  
    { 
#ifdef G4VIS_USE 
      G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();
#endif

#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
      ui->SessionStart();
      delete ui;
#endif

#ifdef G4VIS_USE
     delete visManager;
#endif
    }
  else
    // Batch mode
    {     
      G4String command =("/control/execute ");
      G4String fileName = argv[1];
      G4cout <<"macro --> "<< fileName << G4endl;
      UI->ApplyCommand(command+fileName);
    }

#ifdef G4ANALYSIS_USE
  analysis->finish();
#endif

  delete runManager;

  return 0;
}
