#ifndef StatAccepTestEventAction_h
#define StatAccepTestEventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"

class StatAccepTestSteppingAction;


class StatAccepTestEventAction: public G4UserEventAction {

public:

  StatAccepTestEventAction();
  ~StatAccepTestEventAction();

  virtual void StartOfEventAction(const G4Event* evt);    
  virtual void EndOfEventAction(const G4Event* evt);    

private:

  void instanciateSteppingAction();      

  StatAccepTestSteppingAction* theSteppingAction;

};

#endif
