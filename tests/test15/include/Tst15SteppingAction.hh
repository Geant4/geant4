
#ifndef Tst15SteppingAction_H
#define Tst15SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class Tst15SteppingAction : public G4UserSteppingAction
{
  public:
    Tst15SteppingAction();
    ~Tst15SteppingAction();
    virtual void UserSteppingAction(const G4Step*);
};

#endif

