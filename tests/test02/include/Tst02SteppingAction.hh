
#ifndef Tst02SteppingAction_H
#define Tst02SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class Tst02SteppingAction : public G4UserSteppingAction
{
  public:
    Tst02SteppingAction();
    ~Tst02SteppingAction();
    virtual void UserSteppingAction();
};

#endif

