///////////////////////////////////////////////////////////////////////////////
// File: CCalEndOfEventAction.hh
// Description: CCalEndOfEventAction provides User actions at end of event
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalEndOfEventAction_h
#define CCalEndOfEventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"

class CCaloOrganization;
class G4HCofThisEvent;
class CCalSteppingAction;
class CCalPrimaryGeneratorAction;

typedef G4String nameType;

class CCalEndOfEventAction: public G4UserEventAction{

public:
  CCalEndOfEventAction(CCalPrimaryGeneratorAction*);
  ~CCalEndOfEventAction();
  virtual void StartOfEventAction(const G4Event* evt);    
  virtual void EndOfEventAction(const G4Event* evt);    

private:
  void instanciateSteppingAction();
  void initialize();    

private:
  G4bool                      isInitialized;

  CCalPrimaryGeneratorAction* primaryGenerator;
  CCalSteppingAction*         theSteppingAction;
  nameType*                   SDnames;
  G4int                       numberOfSD;
  CCaloOrganization*          theOrg;

};

#endif
