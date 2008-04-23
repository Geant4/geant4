#ifndef MyEventAction_h
#define MyEventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"

class MyRunAction;
class MySteppingAction;


class MyEventAction: public G4UserEventAction {

public:

  MyEventAction();
  ~MyEventAction();

  virtual void BeginOfEventAction( const G4Event* evt );    
  virtual void EndOfEventAction( const G4Event* evt );    

private:

  void instanciateRunAction();
  void instanciateSteppingAction();

  MyRunAction* theRunAction;
  MySteppingAction* theSteppingAction;

};

#endif
