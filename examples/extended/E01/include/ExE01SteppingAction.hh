// $Id: ExE01SteppingAction.hh,v 1.1 1998/10/14 15:19:40 allison Exp $

#ifndef ExE01SteppingAction_h
#define ExE01SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "Randomize.hh"

#include "histograms.h"

class ExE01SteppingAction : public G4UserSteppingAction {
  private:
    HepJamesRandom theJamesEngine;
    DRand48Engine theDRand48Engine;

  public:
    ExE01SteppingAction(){};
    ~ExE01SteppingAction(){};

    void UserSteppingAction();
};

#endif
