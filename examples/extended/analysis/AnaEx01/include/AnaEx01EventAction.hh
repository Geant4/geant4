// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01EventAction.hh,v 1.2 2000-10-31 13:10:00 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef AnaEx01EventAction_h
#define AnaEx01EventAction_h

#include "G4UserEventAction.hh"

class AnaEx01AnalysisManager;

class AnaEx01EventAction : public G4UserEventAction {
public:
  AnaEx01EventAction(AnaEx01AnalysisManager*);
  virtual ~AnaEx01EventAction();
public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);
private:
  AnaEx01AnalysisManager* fAnalysisManager;
};

#endif

    
