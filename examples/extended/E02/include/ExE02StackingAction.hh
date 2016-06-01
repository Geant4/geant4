
#ifndef ExE02StackingAction_H
#define ExE02StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"
#include <rw/tpordvec.h>

class G4Track;

class ExE02StackingAction : public G4UserStackingAction
{
  public:
    ExE02StackingAction();
    ~ExE02StackingAction();

  public:
    G4ClassificationOfNewTrack ClassifyNewTrack(G4Track *const aTrack);
    void NewStage();
    void PrepareNewEvent();

};

#endif

