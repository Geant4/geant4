
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
    virtual ~ExE02StackingAction();

  public:
   virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
   virtual void NewStage();
   virtual void PrepareNewEvent();

};

#endif

