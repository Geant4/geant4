#ifndef MyRunAction_h
#define MyRunAction_h

#include "G4UserRunAction.hh"

class G4Run;


class MyRunAction: public G4UserRunAction {

public:

  MyRunAction();
  virtual ~MyRunAction();
      
  virtual void BeginOfRunAction( const G4Run* aRun );    
  virtual void EndOfRunAction( const G4Run* aRun );    

  void updateEndOfEvent( const double eventTotalDepositedEnergy );
  // This method is called by MyEventAction, at the end of
  // the event, to provide the total deposited energy in
  // the whole experimental hall.

private:

  int theNumberOfEvents;
  double theSumOfTotalDepositedEnergy;

};

#endif
