#ifndef Tst68EventAction_h
#define Tst68EventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"

class Tst68SteppingAction;
class G4Timer;


class Tst68EventAction: public G4UserEventAction {

public:

  Tst68EventAction();
  ~Tst68EventAction();

  virtual void BeginOfEventAction( const G4Event* evt );    
  virtual void EndOfEventAction( const G4Event* evt );    

private:

  void instanciateSteppingAction();      

  Tst68SteppingAction* theSteppingAction;

  G4Timer* eventTimer;

  G4int numberOfEvents;
  G4double sumTotalDepositedEnergy; 

};

#endif
