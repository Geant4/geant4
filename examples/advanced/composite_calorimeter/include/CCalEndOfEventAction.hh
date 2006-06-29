//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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
