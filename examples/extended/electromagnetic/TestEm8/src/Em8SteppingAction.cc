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
//
// $Id: Em8SteppingAction.cc,v 1.10 2007/11/12 10:54:49 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em8DetectorConstruction.hh"
//#include "G4EnergyLossTables.hh"
//#include "G4SteppingManager.hh"
//#include "G4TrackVector.hh"
#include "Em8SteppingAction.hh"
#include "Em8PrimaryGeneratorAction.hh"
#include "Em8EventAction.hh"
#include "Em8RunAction.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
//#include "G4EventManager.hh"
//#include "G4ios.hh"
//#include <iomanip>
//#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em8SteppingAction::Em8SteppingAction(Em8DetectorConstruction* DET,
                                     Em8EventAction* EA,
                                     Em8RunAction* RA)
:detector (DET),eventaction (EA),runaction (RA),
 IDold(-1) ,evnoold(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em8SteppingAction::~Em8SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em8SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
  G4Track* track = aStep->GetTrack();
  G4int trackID = track->GetTrackID();
  G4int stepID = track->GetCurrentStepNumber();
  if(trackID == 1 && stepID == 1)
    runaction->Fillvertexz(track->GetVertexPosition().z());

  const G4ParticleDefinition* part = track->GetDynamicParticle()->GetDefinition();
  G4VPhysicalVolume* preVolume = aStep->GetPreStepPoint()->GetPhysicalVolume();  
  G4VPhysicalVolume* postVolume = aStep->GetPostStepPoint()->GetPhysicalVolume();  
  G4double ekin = track->GetKineticEnergy() ;  

  // Secondary - first step
  if(trackID != 1 && stepID == 1) {
     if (part == G4Gamma::Gamma() ) {
       eventaction->AddNeutral();

     } else if (part == G4Electron::Electron()) {
       eventaction->AddE();
       runaction->FillTsec(ekin);

     } else if (part == G4Positron::Positron()) {
       eventaction->AddP();
       runaction->FillTsec(ekin);
     }
  }

  // step inside absorber
  if(preVolume == postVolume) { 
    if(preVolume->GetName() == "Absorber") {
      if(part->GetPDGCharge() == 0.0) eventaction->CountStepsNeutral();
      else eventaction->CountStepsCharged();
    }
    // exit of absorber
  } else if (preVolume->GetName() == "Absorber") {

    // primary
    if (trackID == 1) {
      if (track->GetMomentumDirection().z()>0.) {
	eventaction->SetTr();
	runaction->FillTh(track->GetMomentumDirection().theta()) ;
	runaction->FillTt(ekin) ;
	G4double xend= track->GetPosition().x() ;
	G4double yend= track->GetPosition().y() ;
	runaction->FillR(std::sqrt(xend*xend+yend*yend));
      } else {
	eventaction->SetRef();
	runaction->FillThBack(track->GetMomentumDirection().theta() - 0.5*pi) ;
	runaction->FillTb(ekin);
      }
    }
    // transmitted gamma
    if(track->GetMomentumDirection().z()>0. && part == G4Gamma::Gamma()) {
      runaction->FillGammaSpectrum(ekin);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

