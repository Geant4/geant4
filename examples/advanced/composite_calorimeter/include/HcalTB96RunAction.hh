#ifndef HcalTB96RunAction_h
#define HcalTB96RunAction_h
//
// Veronique Lefebure
// 15.11.99
//
#include "G4UserRunAction.hh"

class G4Run;
class HcalTB96RunAction: public G4UserRunAction{

public:
  HcalTB96RunAction(){};
  virtual ~HcalTB96RunAction(){};
      
public:
  virtual void StartOfRunAction(const G4Run* aRun);    
  virtual void EndOfRunAction(const G4Run* aRun);    


};

#endif
