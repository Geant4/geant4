// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01SteppingAction.hh,v 1.1.1.1 2000-09-14 11:37:21 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef AnaEx01SteppingAction_h
#define AnaEx01SteppingAction_h

#include "G4UserSteppingAction.hh"

class G4AnalysisManager;

class AnaEx01SteppingAction : public G4UserSteppingAction {
public:
  AnaEx01SteppingAction(G4AnalysisManager*);
  virtual ~AnaEx01SteppingAction();
  virtual void UserSteppingAction(const G4Step*);
private:
  G4AnalysisManager* fAnalysisManager;
};

#endif
