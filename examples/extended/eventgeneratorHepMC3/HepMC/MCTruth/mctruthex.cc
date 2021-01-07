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
/// \file eventgenerator/HepMC/MCTruth/mctruthex.cc
/// \brief Main program of the eventgenerator/HepMC/MCTruth example
//
//
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "FTFP_BERT.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "MCTruthTrackingAction.hh"
#include "MCTruthEventAction.hh"

#include "MCTruthManager.hh"

int main()
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction);
  runManager->SetUserInitialization(new FTFP_BERT);

  // set mandatory user action class
  runManager->SetUserAction(new PrimaryGeneratorAction);

  // set MCTruth user action classes
  runManager->SetUserAction(new MCTruthTrackingAction);
  runManager->SetUserAction(new MCTruthEventAction);

  // Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/run/verbose 1");
  UI->ApplyCommand("/event/verbose 1");
  UI->ApplyCommand("/tracking/verbose 1");

  // configure MCTruth handling
  MCTruthConfig* config = new MCTruthConfig;
  config->SetMinE(1000.0);
  config->AddParticleType(11);
  MCTruthManager::GetInstance()->SetConfig(config);

  // start a run
  int numberOfEvent = 1;
  runManager->BeamOn(numberOfEvent);

  // job termination
  delete runManager;
}


