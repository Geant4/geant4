// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01EventAction.cc,v 1.1.1.1 2000-09-14 11:37:21 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "G4AnalysisManager.hh"

#include "AnaEx01EventAction.hh"

AnaEx01EventAction::AnaEx01EventAction(
 G4AnalysisManager* aAnalysisManager
):fAnalysisManager(aAnalysisManager){}

AnaEx01EventAction::~AnaEx01EventAction(){}

void AnaEx01EventAction::BeginOfEventAction(const G4Event* aEvent){
  if(fAnalysisManager) fAnalysisManager->BeginOfEvent(aEvent);
}

void AnaEx01EventAction::EndOfEventAction(const G4Event* aEvent) {
  if(fAnalysisManager) fAnalysisManager->EndOfEvent(aEvent);
}

