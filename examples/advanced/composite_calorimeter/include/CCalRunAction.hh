///////////////////////////////////////////////////////////////////////////////
// File: CCalRunAction.hh
// Description: A class for providing user actions at begin and end of run
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalRunAction_h
#define CCalRunAction_h

#include "G4UserRunAction.hh"

class G4Run;
class CCalRunAction: public G4UserRunAction{

public:
  CCalRunAction(){};
  virtual ~CCalRunAction(){};
      
public:
  virtual void BeginOfRunAction(const G4Run* aRun);    
  virtual void EndOfRunAction(const G4Run* aRun);    


};

#endif
