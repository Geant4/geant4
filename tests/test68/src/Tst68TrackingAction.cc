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
#include "Tst68TrackingAction.hh"
#include "G4Track.hh"


Tst68TrackingAction::Tst68TrackingAction() {}


Tst68TrackingAction::~Tst68TrackingAction() {}


void Tst68TrackingAction::PreUserTrackingAction( const G4Track* aTrack ) {
  // This method is called not only once when a particle is created, 
  // but also each time it is resumed, in the case the track gets
  // suspended, as it happens in the case of neutrons with _HP
  // Physics Lists. To be sure that we collect information about
  // a track one and only once, we require that the current step
  // be the first one.
  if ( aTrack->GetCurrentStepNumber() == 0 ) {
    //G4cout << " TrackingAction: id=" << aTrack->GetTrackID()
    //       << "\t" << aTrack->GetDefinition()->GetParticleName()
    //       << "\t" << aTrack->GetKineticEnergy() << " MeV" << G4endl;
  } else {
    //G4cout << " aTrack->GetCurrentStepNumber() = " 
    //	     <<  aTrack->GetCurrentStepNumber() << G4endl;
  }
}


void Tst68TrackingAction::PostUserTrackingAction( const G4Track* aTrack ) {
  // This method is called not only once at the end of the life of
  // a track, but also each time it is suspended, as it happens 
  // in the case of neutrons with _HP Physics Lists. 
  // To be sure that we collect information about a track only once
  // when its life comes to the end, we have to require that its 
  // status is "fStopAndKill".
  if ( aTrack->GetTrackStatus() == fStopAndKill ) {
    //G4cout << " POST-USER Tracking Action info : " << G4endl
    //       << "\t track id = " << aTrack->GetTrackID() << G4endl
    //	     << "\t particle = " << aTrack->GetDefinition()->GetParticleName() << G4endl
    //	     << "\t track status = " << aTrack->GetTrackStatus() << G4endl
    //	     << "\t track length = " << aTrack->GetTrackLength() / mm << " mm" << G4endl
    //	     << "\t track Ekin = " << aTrack->GetKineticEnergy() << " MeV" << G4endl
    //	     << "\t volume = " << aTrack->GetVolume()->GetName() << G4endl
    //	     << "\t material = " << aTrack->GetMaterial()->GetName() << G4endl;
  }
}

