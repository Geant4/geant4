
#ifndef Tst14SteppingAction_H
#define Tst14SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class Tst14SteppingAction : public G4UserSteppingAction
{
  public:
    Tst14SteppingAction();
    virtual ~Tst14SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
};

#endif

