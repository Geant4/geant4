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
// $Id: EventAction.cc,v 1.1 2003-07-31 01:21:16 dwright Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "EventAction.hh"

#include "RunAction.hh"
#include "EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"


EventAction::EventAction(RunAction* runAction)
  :printModulo(1), theRunAction(runAction)
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

  particleList.clear();
}


void EventAction::EndOfEventAction(const G4Event* evt)
{
  // Extract info from stored particles

  G4double kineticEnergy = 0;
  G4double charge = 0;
  G4double cosTheta = 0;

  for (G4int i = 0; i < (G4int)particleList.size(); i++) {
    if ( particleList[i]->GetDefinition()->GetParticleName() == "pi+"  ||
         particleList[i]->GetDefinition()->GetParticleName() == "pi-" )
      {
        kineticEnergy = particleList[i]->GetKineticEnergy();
        charge = particleList[i]->GetCharge();
        cosTheta = particleList[i]->GetMomentumDirection().z();
        theRunAction->FillHists(charge, kineticEnergy, cosTheta);
      }
  }
} 


void EventAction::StoreDynamicParticle(const G4DynamicParticle* dp)
{
  // Copy particle, then store it
 
  G4DynamicParticle* cdp = new G4DynamicParticle(*dp);
  particleList.push_back(cdp);
}







