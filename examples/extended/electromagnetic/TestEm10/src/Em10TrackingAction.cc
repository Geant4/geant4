 #include "Em10TrackingAction.hh"

 #include "G4TrackingManager.hh"
 #include "G4Track.hh"
 #include "G4RunManager.hh"

 #include "G4UImanager.hh"

 Em10TrackingAction::Em10TrackingAction() 
 { }

 void Em10TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
 {

  if( aTrack->GetParentID() == 1 && aTrack->GetKineticEnergy() > 100.*GeV ){
    G4cout << "[Em10TrackingAction::DEBUG]" << G4endl;
    G4cout << " Track ID:          " << aTrack->GetTrackID() << G4endl;
    G4cout << " particle:          " << aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName() << G4endl;
    G4cout << " Parent ID:         " << aTrack->GetParentID() << G4endl;
    G4cout << " created by:        " << aTrack->GetCreatorProcess()->GetProcessName() << G4endl;
    G4cout << " kin. energy (TeV): " << aTrack->GetKineticEnergy() / TeV << G4endl;
    G4cout << " volume:            " << aTrack->GetVolume()->GetName() << G4endl;
    G4cout << " global time:       " << aTrack->GetGlobalTime() << G4endl;

    G4cout << " Killing event..." << G4endl;
    if( aTrack->GetTrackID() != 1 )
      const_cast<G4Track*>(aTrack)->SetTrackStatus( fKillTrackAndSecondaries );
  }

  return;
 }
