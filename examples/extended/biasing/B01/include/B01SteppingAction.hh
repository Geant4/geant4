#ifndef B01SteppingAction_hh
#define B01SteppingAction_hh B01SteppingAction_hh

#include "G4UserSteppingAction.hh"

class B01SteppingAction : public G4UserSteppingAction {
public:
  B01SteppingAction();
  virtual ~B01SteppingAction();
  virtual void UserSteppingAction(const G4Step* pStep);
};

  
#endif
