#include "B01SteppingAction.hh"
#include "G4Track.hh"

B01SteppingAction::B01SteppingAction()
{}

B01SteppingAction::~B01SteppingAction()
{}


void B01SteppingAction::UserSteppingAction(const G4Step* pStep) {
  G4Track* pTrack = pStep->GetTrack();
  G4ParticleDefinition* pParticle = pTrack->GetDefinition();
  G4String szName = pParticle->GetParticleName();
  if(szName=="alpha" || 
     szName=="deuteron" ||
     szName=="proton" ||
     szName=="neutron" || 
     szName=="gamma"){
    return;
  }
  pTrack->SetTrackStatus(fStopAndKill);
}
  
