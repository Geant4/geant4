// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01SteppingAction.hh,v 1.3 2000-11-10 14:12:07 gbarrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef AnaEx01SteppingAction_h
#define AnaEx01SteppingAction_h

#include "G4UserSteppingAction.hh"

class AnaEx01AnalysisManager;

class AnaEx01SteppingAction : public G4UserSteppingAction {
public:
  AnaEx01SteppingAction(AnaEx01AnalysisManager* = 0);
  virtual ~AnaEx01SteppingAction();
  virtual void UserSteppingAction(const G4Step*);
private:
  AnaEx01AnalysisManager* fAnalysisManager;
};

#endif
