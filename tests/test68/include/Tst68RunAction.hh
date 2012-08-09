#ifndef Tst68RunAction_h
#define Tst68RunAction_h

#include "G4UserRunAction.hh"

class G4Run;


class Tst68RunAction: public G4UserRunAction {

public:

  Tst68RunAction(){};
  virtual ~Tst68RunAction();
      
  virtual void BeginOfRunAction( const G4Run* aRun );    
  virtual void EndOfRunAction( const G4Run* aRun );    

};

#endif
