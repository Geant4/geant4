
#ifndef MyStackingAction_H
#define MyStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"
#include <rw/tpordvec.h>

class G4Track;

class MyStackingAction : public G4UserStackingAction
{
  public:
    MyStackingAction();
    ~MyStackingAction();

  public:
    G4ClassificationOfNewTrack ClassifyNewTrack(G4Track *const aTrack);
    void NewStage();
    void PrepareNewEvent();

};

#endif

