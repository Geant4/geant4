// $Id: ExE01SteppingAction.hh,v 1.2 1999/04/17 04:06:43 kurasige Exp $

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
    virtual ~ExE01SteppingAction(){};

    virtual void UserSteppingAction(const G4Step*);
};

#endif
