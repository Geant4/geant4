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
// $Id: test50.cc,v 1.27 2003-06-25 10:20:11 gunter Exp $
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
#include "Tst50PhysicsList.hh"
#include "Tst50PrimaryGeneratorAction.hh"
#include "Tst50RunAction.hh"
#include "Tst50EventAction.hh"
#include "Tst50SteppingAction.hh"
#include "Tst50SteppingVerbose.hh"
#include "Tst50AnalysisManager.hh"
#ifdef G4VIS_USE
#include "Tst50VisManager.hh"
#endif

int main(int argc,char** argv) {
 	
  HepRandom::setTheEngine(new RanecuEngine);   
  G4int seed = time(0);
  HepRandom ::setTheSeed(seed);
 
  G4RunManager * runManager = new G4RunManager;

  Tst50DetectorConstruction* tst50Detector = new Tst50DetectorConstruction();
  runManager->SetUserInitialization(tst50Detector);  

  Tst50PhysicsList* tst50Physics = new Tst50PhysicsList();
  runManager->SetUserInitialization(tst50Physics);

#ifdef G4VIS_USE 
  G4VisManager* visManager = new Tst50VisManager;
  visManager->Initialize();
#endif

  Tst50PrimaryGeneratorAction* tst50PrimaryParticle = new Tst50PrimaryGeneratorAction(); 
  runManager->SetUserAction(tst50PrimaryParticle);

  Tst50RunAction* tst50Run=new Tst50RunAction(); 
  runManager->SetUserAction(tst50Run);  

  Tst50EventAction *tst50EventAction = new Tst50EventAction();
 
  runManager->SetUserAction(tst50EventAction);
     
  Tst50SteppingAction* tst50SteppingAction = new Tst50SteppingAction(tst50PrimaryParticle,
                                                                tst50Run,
                                                                tst50Detector);
  runManager->SetUserAction(tst50SteppingAction);

  Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
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
