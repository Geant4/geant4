#ifndef ExN06StackingAction_H
#define ExN06StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class ExN06StackingAction : public G4UserStackingAction
{
public:
  ExN06StackingAction();
  virtual ~ExN06StackingAction();

public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();

private:
  G4int gammaCounter;
};

#endif

