
#ifndef Tst13SteppingAction_H
#define Tst13SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class Tst13SteppingAction : public G4UserSteppingAction
{
  public:
    Tst13SteppingAction();
    ~Tst13SteppingAction();
    virtual void UserSteppingAction();
};

#endif

