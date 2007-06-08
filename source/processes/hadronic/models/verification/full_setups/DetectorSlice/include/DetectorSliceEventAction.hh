#ifndef DetectorSliceEventAction_h
#define DetectorSliceEventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"

class DetectorSliceSteppingAction;
class G4Timer;


class DetectorSliceEventAction: public G4UserEventAction {

public:

  DetectorSliceEventAction();
  ~DetectorSliceEventAction();

  virtual void BeginOfEventAction( const G4Event* evt );    
  virtual void EndOfEventAction( const G4Event* evt );    

private:

  void instanciateSteppingAction();      

  DetectorSliceSteppingAction* theSteppingAction;

  G4Timer* eventTimer;

};

#endif
