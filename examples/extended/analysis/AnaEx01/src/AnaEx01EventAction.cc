// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01EventAction.cc,v 1.3 2000-11-10 14:08:08 gbarrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifdef G4ANALYSIS_USE
#include "AnaEx01AnalysisManager.hh"
#endif

#include "AnaEx01EventAction.hh"

AnaEx01EventAction::AnaEx01EventAction(
 AnaEx01AnalysisManager* aAnalysisManager
):fAnalysisManager(aAnalysisManager){}

AnaEx01EventAction::~AnaEx01EventAction(){}

void AnaEx01EventAction::BeginOfEventAction(const G4Event* aEvent){
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->BeginOfEvent(aEvent);
#endif
}

void AnaEx01EventAction::EndOfEventAction(const G4Event* aEvent) {
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->EndOfEvent(aEvent);
#endif
}

