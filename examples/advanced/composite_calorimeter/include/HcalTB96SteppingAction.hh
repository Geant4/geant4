#ifndef HcalTB96SteppingAction_h
#define HcalTB96SteppingAction_h
//
// Veronique Lefebure
// 15.11.99
//

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4Step;
class HcalTB96SteppingAction : public G4UserSteppingAction {

  friend class HcalTB96EndOfEventAction;
private:
  HcalTB96SteppingAction();
public:
  ~HcalTB96SteppingAction();
      
public:
  virtual void UserSteppingAction(const G4Step* aStep);    


private:

  float          LateralProfile[28];	
  float          timeDeposit[50];
  int            timeHistoMaxBin;	    
      
private:
  void endOfEvent();
};

#endif
