#ifndef MyRunAction_h
#define MyRunAction_h

#include "G4UserRunAction.hh"

class G4Run;


class MyRunAction: public G4UserRunAction {

public:

  MyRunAction(){};
  virtual ~MyRunAction(){};
      
  virtual void BeginOfRunAction(const G4Run* aRun);    
  virtual void EndOfRunAction(const G4Run* aRun);    

};

#endif
