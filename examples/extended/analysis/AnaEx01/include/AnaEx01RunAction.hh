// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01RunAction.hh,v 1.1.1.1 2000-09-14 11:37:21 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef AnaEx01RunAction_h
#define AnaEx01RunAction_h

#include "G4UserRunAction.hh"

class G4AnalysisManager;

class AnaEx01RunAction : public G4UserRunAction {
public:
  AnaEx01RunAction(G4AnalysisManager*);
  ~AnaEx01RunAction();
public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);  
private:
  G4AnalysisManager* fAnalysisManager;
};

#endif

