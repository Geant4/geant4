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
///////////////////////////////////////////////////////////////////////////////
// File: CCalRunAction.cc
// Description: A class for providing user actions at begin and end of run
///////////////////////////////////////////////////////////////////////////////
#include "CCalRunAction.hh"

#include "globals.hh"
#include "G4Run.hh"

#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

#include "CCalAnalysis.hh"


void CCalRunAction::BeginOfRunAction(const G4Run* aRun) {

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // A.R. Added for visualization of events.
  if ( G4VVisManager::GetConcreteInstance() ) {
    G4UImanager* UI = G4UImanager::GetUIpointer(); 
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  } 

  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  analysis->BeginOfRun(aRun->GetRunID());

  
}


void CCalRunAction::EndOfRunAction(const G4Run* aRun) {

  G4cout << "### Run " << aRun->GetRunID() << " end." << G4endl;

  // A.R. Added for visualization of events.
  if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  analysis->EndOfRun(aRun->GetRunID());

}   

