#ifndef RandomCaloRunAction_h
#define RandomCaloRunAction_h

#include "G4UserRunAction.hh"

class G4Run;


class RandomCaloRunAction: public G4UserRunAction {

public:

  RandomCaloRunAction() {};
  virtual ~RandomCaloRunAction();
      
  virtual void BeginOfRunAction( const G4Run* aRun );    
  virtual void EndOfRunAction( const G4Run* aRun );    

};

#endif
