
#ifndef SteppingAction_H
#define SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction();
    virtual ~SteppingAction();

    virtual void UserSteppingAction(const G4Step* );
};

#endif

