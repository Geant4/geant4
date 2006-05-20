#ifndef StatAccepTestRunAction_h
#define StatAccepTestRunAction_h

#include "G4UserRunAction.hh"

class G4Run;


class StatAccepTestRunAction: public G4UserRunAction {

public:

  StatAccepTestRunAction(){};
  virtual ~StatAccepTestRunAction();
      
  virtual void BeginOfRunAction( const G4Run* aRun );    
  virtual void EndOfRunAction( const G4Run* aRun );    

};

#endif
