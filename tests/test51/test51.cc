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
// $Id: test51.cc,v 1.1 2005-07-05 11:09:44 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "globals.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Tst51DetectorConstruction.hh"
#include "Tst51PhysicsList.hh"
#include "Tst51PrimaryGeneratorAction.hh"
#include "Tst51RunAction.hh"
#include "Tst51EventAction.hh"
#include "Tst51SteppingAction.hh"
#include "Tst51SteppingVerbose.hh"
#include "Tst51AnalysisManager.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

int main(int argc,char** argv) {
 	
  HepRandom :: setTheEngine(new RanecuEngine);   
  G4int seed = time(0);
  HepRandom :: setTheSeed(seed);
 
  G4RunManager * runManager = new G4RunManager;

  Tst51DetectorConstruction* tst51Detector = new Tst51DetectorConstruction();
  runManager->SetUserInitialization(tst51Detector);  
  G4cout<<"Experimental setup loaded"<<G4endl;
  
  Tst51PhysicsList* tst51Physics = new Tst51PhysicsList();
  runManager->SetUserInitialization(tst51Physics);
  G4cout<<"Physics loaded"<<G4endl;

#ifdef G4VIS_USE 
  G4VisExecutive* visManager = new G4VisExecutive;
  visManager -> Initialize();  
  G4cout<<"Visualisation loaded"<<G4endl;
#endif

  Tst51PrimaryGeneratorAction* tst51PrimaryParticle = new Tst51PrimaryGeneratorAction(); 
  runManager->SetUserAction(tst51PrimaryParticle);
  G4cout<<"Primary particle loaded"<<G4endl;

  Tst51RunAction* tst51Run=new Tst51RunAction(); 
  runManager->SetUserAction(tst51Run);  
  G4cout<<"RunAction loaded"<<G4endl;

  Tst51EventAction *tst51EventAction = new Tst51EventAction();
 
  runManager->SetUserAction(tst51EventAction);
  G4cout<<"EventAction loaded"<<G4endl;  
  
  Tst51SteppingAction* tst51SteppingAction = new Tst51SteppingAction();
  runManager->SetUserAction(tst51SteppingAction);
  G4cout<<"SteppingAction loaded"<<G4endl;

  Tst51AnalysisManager* analysis = Tst51AnalysisManager::getInstance();
  analysis->book();
  G4cout<<"Analysis loaded"<<G4endl;

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

      UI->ApplyCommand("/control/execute run.mac");
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
