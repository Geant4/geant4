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
#define hTestTrackingAction_CPP

//---------------------------------------------------------------------------
//
// ClassName:   hTestTrackingAction
//  
// Description: Implementation file for MC truth.
//
// Author:      V.Ivanchenko 17/03/01
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestTrackingAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestTrackingAction::hTestTrackingAction():
  theHisto(hTestHisto::GetPointer())
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestTrackingAction::~hTestTrackingAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

  G4ParticleDefinition* particle = aTrack->GetDefinition();
  G4String name = particle->GetParticleName();
  theHisto->ResetTrackLength();
  G4int pid = aTrack->GetParentID(); 

  if(0 < theHisto->GetVerbose() &&
    (theHisto->GetMaxEnergy() < aTrack->GetKineticEnergy() && pid > 0)) {
    G4cout << "Track #" 
           << aTrack->GetTrackID() << " of " << name
           << " Emax(MeV)= " << theHisto->GetMaxEnergy()/MeV
           << " Ekin(MeV)= " << aTrack->GetKineticEnergy()/MeV
           << " ## EventID= " 
           << (G4EventManager::GetEventManager())->GetConstCurrentEvent()
              ->GetEventID() 
           << G4endl;
  }

  //Save primary parameters

  if(0 == pid) {

    //    theHisto->StartEvent();      

    G4double kinE = aTrack->GetKineticEnergy();
    theHisto->SaveToTuple("TKIN", kinE/MeV);      

    G4double mass = 0.0;
    if(particle) {
      mass = particle->GetPDGMass();
      theHisto->SaveToTuple("MASS", mass/MeV);      
      theHisto->SaveToTuple("CHAR",(particle->GetPDGCharge())/eplus);      
      G4double beta = 1.;
	if(mass > 0.) {
          G4double gamma = kinE/mass + 1.;
          beta = sqrt(1. - 1./(gamma*gamma));
	}
      theHisto->SaveToTuple("BETA", beta);            
    }

    G4ThreeVector pos = aTrack->GetPosition();
    G4ThreeVector dir = aTrack->GetMomentumDirection();
    /*
    theHisto->SaveToTuple("X0",(pos.x())/mm);
    theHisto->SaveToTuple("Y0",(pos.y())/mm);
    theHisto->SaveToTuple("Z0",(pos.z())/mm);
    */
    if(1 < theHisto->GetVerbose()) {
      G4cout << "hTestTrackingAction: Primary kinE(MeV)= " << kinE/MeV
           << "; m(MeV)= " << mass/MeV   
           << "; pos= " << pos << ";  dir= " << dir << G4endl;
    }

    // delta-electron
  } else if ("e-" == name) {
    if(1 < theHisto->GetVerbose()) {
      G4cout << "hTestTrackingAction: Secondary electron " << G4endl;
    }
    theHisto->AddDeltaElectron(aTrack->GetDynamicParticle());

  } else if ("e+" == name) {
    if(1 < theHisto->GetVerbose()) {
      G4cout << "hTestTrackingAction: Secondary positron " << G4endl;
    }
    //    theHisto->AddPositron(aTrack->GetDynamicParticle());

  } else if ("gamma" == name) {
    if(0 < theHisto->GetVerbose()) {
      G4cout << "hTestTrackingAction: Secondary gamma " 
             << aTrack->GetDynamicParticle() << G4endl;
      aTrack->GetDynamicParticle()->DumpInfo();
    }
    theHisto->AddPhoton(aTrack->GetDynamicParticle());
 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....








