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
//
//
// $Id: EventAction.cc,v 1.1 2003-05-27 13:44:48 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "EventAction.hh"

#include "TargetHit.hh"
#include "RunAction.hh"
#include "EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"


EventAction::EventAction(RunAction* runAction)
  :collID(-1), printModulo(1), theRunAction(runAction)
{
  eventMessenger = new EventActionMessenger(this);
}


EventAction::~EventAction()
{
  delete eventMessenger;
}


void EventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();
 
  if (evtNb%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    HepRandom::showEngineStatus();
  }

  if (collID==-1) {
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    collID = SDman->GetCollectionID("TgtCollection");
  }
}


void EventAction::EndOfEventAction(const G4Event* evt)
{
  // Extract info from hits

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  TargetHitsCollection* THC = 0;

  G4int n_hit = 0;
  G4double cosTheta = 0;
  G4double kineticEnergy = 0;
  G4double charge = 0;

  if (HCE) THC = (TargetHitsCollection*)(HCE->GetHC(collID));

  if (THC) {
    n_hit = THC->entries();
    for (G4int i=0;i<n_hit;i++) {
      cosTheta = (*THC)[i]->GetTheta(); 
      kineticEnergy = (*THC)[i]->GetEkin();
      charge = (*THC)[i]->GetCharge();
      //      G4cout << "Charge, KE, cos(theta)" << charge << " , "
      //             << kineticEnergy << " , " << cosTheta << G4endl;

      theRunAction->FillHists(charge, kineticEnergy, cosTheta);
    }
  }
} 









