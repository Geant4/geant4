
#ifndef PersEx01StackingAction_H
#define PersEx01StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"

class G4Track;

class PersEx01StackingAction : public G4UserStackingAction
{
  public:
    PersEx01StackingAction();
    virtual ~PersEx01StackingAction();

  public:
   virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
   virtual void NewStage();
   virtual void PrepareNewEvent();

};

#endif

