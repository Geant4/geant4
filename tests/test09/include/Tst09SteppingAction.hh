
#ifndef Tst09SteppingAction_H
#define Tst09SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class Tst09SteppingAction : public G4UserSteppingAction
{
  public:
    Tst09SteppingAction();
    ~Tst09SteppingAction();
    virtual void UserSteppingAction();
};

#endif

