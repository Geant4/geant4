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

#ifdef G4ANALYSIS_USE
#include "CCalAnalysis.hh"
#endif


void CCalRunAction::BeginOfRunAction(const G4Run* aRun) {

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // A.R. Added for visualization of events.
  if ( G4VVisManager::GetConcreteInstance() ) {
    G4UImanager* UI = G4UImanager::GetUIpointer(); 
    UI->ApplyCommand("/vis/scene/notifyHandlers");
  } 

#ifdef G4ANALYSIS_USE  
  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  analysis->BeginOfRun(aRun->GetRunID());
#endif
  
}


void CCalRunAction::EndOfRunAction(const G4Run* aRun) {

  G4cout << "Executing CCalRunAction" << G4endl;

  // A.R. Added for visualization of events.
  if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

#ifdef G4ANALYSIS_USE  
  CCalAnalysis* analysis = CCalAnalysis::getInstance();
  analysis->EndOfRun(aRun->GetRunID());
#endif

}   

