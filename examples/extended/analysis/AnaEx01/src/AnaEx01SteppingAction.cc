// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01SteppingAction.cc,v 1.3 2000-11-10 14:08:10 gbarrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifdef G4ANALYSIS_USE
#include "AnaEx01AnalysisManager.hh"
#endif

#include "AnaEx01SteppingAction.hh"

AnaEx01SteppingAction::AnaEx01SteppingAction(
 AnaEx01AnalysisManager* aAnalysisManager
):fAnalysisManager(aAnalysisManager){}

AnaEx01SteppingAction::~AnaEx01SteppingAction(){}
void AnaEx01SteppingAction::UserSteppingAction(const G4Step* aStep){
#ifdef G4ANALYSIS_USE
  if(fAnalysisManager) fAnalysisManager->Step(aStep);
#endif
}



