
#ifndef Tst16SteppingAction_H
#define Tst16SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class Tst16SteppingAction : public G4UserSteppingAction
{
  public:
    Tst16SteppingAction();
    ~Tst16SteppingAction();
    virtual void UserSteppingAction(const G4Step*);
};

#endif

