//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
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
