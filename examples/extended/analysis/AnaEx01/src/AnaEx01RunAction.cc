// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01RunAction.cc,v 1.3 2000-11-10 14:08:09 gbarrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifdef G4ANALYSIS_USE
#include "AnaEx01AnalysisManager.hh"
#endif

#include "AnaEx01RunAction.hh"

AnaEx01RunAction::AnaEx01RunAction(
 AnaEx01AnalysisManager* aAnalysisManager
):fAnalysisManager(aAnalysisManager){}

AnaEx01RunAction::~AnaEx01RunAction(){}

void AnaEx01RunAction::BeginOfRunAction(const G4Run* aRun) {
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->BeginOfRun(aRun);
#endif
}

void AnaEx01RunAction::EndOfRunAction(const G4Run* aRun){
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->EndOfRun(aRun);
#endif
}

