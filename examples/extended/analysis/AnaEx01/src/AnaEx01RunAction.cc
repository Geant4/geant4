// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01RunAction.cc,v 1.2 2000-10-31 13:10:01 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "AnaEx01AnalysisManager.hh"

#include "AnaEx01RunAction.hh"

AnaEx01RunAction::AnaEx01RunAction(
 AnaEx01AnalysisManager* aAnalysisManager
):fAnalysisManager(aAnalysisManager){}

AnaEx01RunAction::~AnaEx01RunAction(){}

void AnaEx01RunAction::BeginOfRunAction(const G4Run* aRun) {
  if(fAnalysisManager) fAnalysisManager->BeginOfRun(aRun);
}

void AnaEx01RunAction::EndOfRunAction(const G4Run* aRun){
  if(fAnalysisManager) fAnalysisManager->EndOfRun(aRun);
}

