
#ifndef MyStackingAction_H
#define MyStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"
#include "g4rw/tpordvec.h"

class G4Track;

class MyStackingAction : public G4UserStackingAction
{
  public:
    MyStackingAction();
    virtual ~MyStackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track *);
    virtual void NewStage();
    virtual void PrepareNewEvent();

};

#endif

