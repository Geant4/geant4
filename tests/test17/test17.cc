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
// $Id: test17.cc,v 1.13 2008-06-24 13:46:00 cirrone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - test17
//
// --------------------------------------------------------------
// Comments
//
// 18.08.2000 V.Ivanchenko clean up visualisation and dummy output
//
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"

#include "Test17DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "Test17PrimaryGeneratorAction.hh"
#include "Test17RunAction.hh"
#include "Test17EventAction.hh"
#include "Test17SteppingAction.hh"

int main(int argc,char** argv) {

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Test17DetectorConstruction* detector;
  detector = new Test17DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList());

  // set user action classes
  runManager->SetUserAction(new Test17PrimaryGeneratorAction(detector));
  Test17RunAction* runaction = new Test17RunAction;
  runManager->SetUserAction(runaction);

  Test17EventAction* eventaction = new Test17EventAction(runaction);
  runManager->SetUserAction(eventaction);

  Test17SteppingAction* steppingaction = new Test17SteppingAction(detector,
                                               eventaction, runaction);
  runManager->SetUserAction(steppingaction);

  //Initialize G4 kernel
  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (argc==1)   // Define UI terminal for interactive mode
    {
      G4UIsession * session = new G4UIterminal;
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  // job termination
  delete runManager;

  return 0;
}

