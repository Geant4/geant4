
#ifndef Tst11SteppingAction_H
#define Tst11SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class Tst11SteppingAction : public G4UserSteppingAction
{
  public:
    Tst11SteppingAction();
    ~Tst11SteppingAction();
    virtual void UserSteppingAction();
};

#endif

