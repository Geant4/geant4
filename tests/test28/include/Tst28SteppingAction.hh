
#ifndef Tst28SteppingAction_H
#define Tst28SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class Tst28SteppingAction : public G4UserSteppingAction
{
  public:
    Tst28SteppingAction();
    virtual ~Tst28SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
};

#endif

