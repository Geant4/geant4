///////////////////////////////////////////////////////////////////////////////
// File: HCalTB96EndOfEventAction.hh
// Date: 15.11.99. 
// Modifications: 09/00 SB
///////////////////////////////////////////////////////////////////////////////

#ifndef HcalTB96EndOfEventAction_h
#define HcalTB96EndOfEventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"

class CCaloOrganization;
class G4HCofThisEvent;
class CCalSteppingAction;
class CMSPrimaryGeneratorAction;

typedef G4String nameType;

class HcalTB96EndOfEventAction: public G4UserEventAction{

public:
  HcalTB96EndOfEventAction(CMSPrimaryGeneratorAction*);
  ~HcalTB96EndOfEventAction();
  virtual void StartOfEventAction(const G4Event* evt);    
  virtual void EndOfEventAction(const G4Event* evt);    

private:
  void instanciateSteppingAction();
  void initialize();    

private:
  G4bool                     isInitialized;

  CMSPrimaryGeneratorAction* primaryGenerator;
  CCalSteppingAction*        theSteppingAction;
  nameType*                  SDnames;
  G4int                      numberOfSD;
  CCaloOrganization*         theOrg;

};

#endif
