
#ifndef Tst12SteppingAction_H
#define Tst12SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class Tst12SteppingAction : public G4UserSteppingAction
{
  public:
    Tst12SteppingAction();
    virtual ~Tst12SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
};

#endif

