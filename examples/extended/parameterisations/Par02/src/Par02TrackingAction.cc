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
/// \file Par02TrackingAction.cc
/// \brief Implementation of the Par02TrackingAction class

#include "Par02TrackingAction.hh"
#include "Par02EventInformation.hh"
#include "Par02PrimaryParticleInformation.hh"
#include "Par02Output.hh"

#include "G4ThreeVector.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackingManager.hh"
#include <iomanip>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02TrackingAction::Par02TrackingAction() : G4UserTrackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02TrackingAction::~Par02TrackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02TrackingAction::PreUserTrackingAction( const G4Track* aTrack ) {
  // Kill the tracks that have a small transverse momentum or that are not
  // in the central region.
  if ( aTrack->GetMomentum().perp() < 1.0*MeV  ||
       std::abs( aTrack->GetMomentum().pseudoRapidity() ) > 5.5 ) {
    ( (G4Track*) aTrack )->SetTrackStatus( fStopAndKill );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02TrackingAction::PostUserTrackingAction( const G4Track* aTrack ) {
  if ( aTrack->GetTrackStatus() == fStopAndKill  &&  aTrack->GetParentID() == 0 ) {
    Par02PrimaryParticleInformation* info = (Par02PrimaryParticleInformation*) 
       aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetUserInformation();
    //info->Print();
    Par02Output::Instance()->SaveTrack( Par02Output::eSaveMC,
                                        info->GetPartID(),
                                        info->GetPDG(),
                                        info->GetMCMomentum()/MeV );
    Par02Output::Instance()->SaveTrack( Par02Output::eSaveTracker,
                                        info->GetPartID(),
                                        info->GetPDG(),
                                        info->GetTrackerMomentum()/MeV,
                                        info->GetTrackerResolution(),
                                        info->GetTrackerEfficiency() );
    Par02Output::Instance()->SaveTrack( Par02Output::eSaveEMCal,
                                        info->GetPartID(),
                                        info->GetPDG(),
                                        info->GetEMCalPosition()/mm,
                                        info->GetEMCalResolution(),
                                        info->GetEMCalEfficiency(),
                                        info->GetEMCalEnergy()/MeV );
    Par02Output::Instance()->SaveTrack( Par02Output::eSaveHCal,
                                        info->GetPartID(),
                                        info->GetPDG(),
                                        info->GetHCalPosition()/mm,
                                        info->GetHCalResolution(),
                                        info->GetHCalEfficiency(),
                                        info->GetHCalEnergy()/MeV );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

