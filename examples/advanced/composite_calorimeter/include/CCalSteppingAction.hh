///////////////////////////////////////////////////////////////////////////////
// File: CCalSeppingAction.hh
// Description: Stepping action needed for analysis
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalSteppingAction_h
#define CCalSteppingAction_h

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4Step;
class CCalSteppingAction : public G4UserSteppingAction {

  friend class CCalEndOfEventAction;
private:
  CCalSteppingAction();
public:
  ~CCalSteppingAction();
      
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
