#ifndef RandomCaloEventAction_h
#define RandomCaloEventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"

class RandomCaloSteppingAction;


class RandomCaloEventAction: public G4UserEventAction {

public:

  RandomCaloEventAction();
  ~RandomCaloEventAction();

  virtual void StartOfEventAction( const G4Event* evt );    
  virtual void EndOfEventAction( const G4Event* evt );    

private:

  void instanciateSteppingAction();      

  RandomCaloSteppingAction* theSteppingAction;

  G4int numberOfEvents;
  G4double sumTotalDepositedEnergy; 

};

#endif
