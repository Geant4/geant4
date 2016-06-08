
#ifndef ExN04SteppingAction_H
#define ExN04SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class ExN04SteppingAction : public G4UserSteppingAction
{
  public:
    ExN04SteppingAction();
    virtual ~ExN04SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
};

#endif

