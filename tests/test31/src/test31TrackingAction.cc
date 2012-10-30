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
#define test31TrackingAction_CPP

//---------------------------------------------------------------------------
//
// ClassName:   test31TrackingAction
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

#include "test31TrackingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31TrackingAction::test31TrackingAction():
  theHisto(test31Histo::GetPointer())
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31TrackingAction::~test31TrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

  G4ParticleDefinition* particle = aTrack->GetDefinition();
  G4String name = particle->GetParticleName();
  theHisto->ResetTrackLength();
  G4int pid = aTrack->GetParentID();

  if(1 < theHisto->GetVerbose() &&
    (theHisto->GetMaxEnergy() < aTrack->GetKineticEnergy() && pid > 0))
  {
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
    theHisto->SaveToTuple("tkin", kinE/MeV);

    G4double mass = 0.0;
    if(particle) {
      mass = particle->GetPDGMass();
      theHisto->SaveToTuple("mass", mass/MeV);
      //theHisto->SaveToTuple("char",(particle->GetPDGCharge())/eplus);
      G4double beta = 1.;
      if(mass > 0.) beta = std::sqrt(kinE*(kinE + 2.0*mass))/(kinE + mass);

      theHisto->SaveToTuple("beta", beta);
    }

    G4ThreeVector pos = aTrack->GetPosition();
    G4ThreeVector dir = aTrack->GetMomentumDirection();
    /*
    theHisto->SaveToTuple("X0",(pos.x())/mm);
    theHisto->SaveToTuple("Y0",(pos.y())/mm);
    theHisto->SaveToTuple("Z0",(pos.z())/mm);
    */
    if(1 < theHisto->GetVerbose()) {
      G4cout << "test31TrackingAction: Primary kinE(MeV)= " << kinE/MeV
           << "; m(MeV)= " << mass/MeV
           << "; pos= " << pos << ";  dir= " << dir << G4endl;
    }

    // delta-electron
  } else if (0 < pid && "e-" == name) {
    if(1 < theHisto->GetVerbose()) {
      G4cout << "test31TrackingAction: Secondary electron " << G4endl;
    }
    theHisto->AddDeltaElectron(aTrack->GetDynamicParticle());

  } else if (0 < pid && "e+" == name) {
    if(1 < theHisto->GetVerbose()) {
      G4cout << "test31TrackingAction: Secondary positron " << G4endl;
    }
    theHisto->AddPositron(aTrack->GetDynamicParticle());

  } else if (0 < pid && "gamma" == name) {
    if(1 < theHisto->GetVerbose()) {
      G4cout << "test31TrackingAction: Secondary gamma " << G4endl;
    }
    theHisto->AddPhoton(aTrack->GetDynamicParticle());
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  if(aTrack->GetParentID() == 0) {
    const G4VProcess* process = 
          aTrack->GetStep()->GetPostStepPoint()->GetProcessDefinedStep();
    if (process) theHisto->CountProcess(process->GetProcessName());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

