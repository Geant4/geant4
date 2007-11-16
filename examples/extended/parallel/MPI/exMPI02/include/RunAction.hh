// $Id: RunAction.hh,v 1.1 2007-11-16 14:29:33 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   RunAction.hh
//
//                                         2007 Q
// ====================================================================
#ifndef RUN_ACTION_H
#define RUN_ACTION_H

#include "G4UserRunAction.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class RunAction : public G4UserRunAction {

public:
  RunAction();
  ~RunAction();

  virtual void BeginOfRunAction(const G4Run* arun);
  virtual void EndOfRunAction(const G4Run* arun);

};

#endif
