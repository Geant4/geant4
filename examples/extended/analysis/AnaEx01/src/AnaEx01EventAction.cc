// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01EventAction.cc,v 1.2 2000-10-31 13:10:01 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "AnaEx01AnalysisManager.hh"

#include "AnaEx01EventAction.hh"

AnaEx01EventAction::AnaEx01EventAction(
 AnaEx01AnalysisManager* aAnalysisManager
):fAnalysisManager(aAnalysisManager){}

AnaEx01EventAction::~AnaEx01EventAction(){}

void AnaEx01EventAction::BeginOfEventAction(const G4Event* aEvent){
  if(fAnalysisManager) fAnalysisManager->BeginOfEvent(aEvent);
}

void AnaEx01EventAction::EndOfEventAction(const G4Event* aEvent) {
  if(fAnalysisManager) fAnalysisManager->EndOfEvent(aEvent);
}

