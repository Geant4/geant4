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
// $Id: test20.cc,v 1.6 2009-08-05 17:33:27 kurasige Exp $
//

#include "Tst20DetectorConstruction.hh"
#include "Tst20RunAction.hh"
#include "Tst20PrimaryGeneratorAction.hh"
#include "Tst20PhysicsList.hh"
#include "Tst20SteppingAction.hh"
#include "Tst20TrackingAction.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  Tst20DetectorConstruction* detector;
  detector = new Tst20DetectorConstruction;

  runManager->SetUserInitialization(detector);
  Tst20PhysicsList* pl;
  runManager->SetUserInitialization(pl = new Tst20PhysicsList);
  pl->DisableCheckParticleList();

  // UserAction classes
  runManager->SetUserAction(new Tst20PrimaryGeneratorAction(detector));
  Tst20RunAction* runaction = new Tst20RunAction;
  runManager->SetUserAction(runaction);

  runManager->SetUserAction(new Tst20SteppingAction);
  //runManager->SetUserAction(new Tst20TrackingAction);

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if(argc==1)
  {
    // G4UIterminal is a (dumb) terminal.
    G4UIsession* session = new G4UIterminal;
    UImanager->ApplyCommand("/control/execute prerunTst20.mac");
    session->SessionStart();
    delete session;
  }
  else
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  delete runManager;
  return 0;
}
