#include "ExN06StackingAction.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"

ExN06StackingAction::ExN06StackingAction()
: gammaCounter(0)
{;}

ExN06StackingAction::~ExN06StackingAction()
{;}

G4ClassificationOfNewTrack
ExN06StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  if(aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition())
  { // particle is optical photon
    if(aTrack->GetParentID()>0)
    { // particle is secondary
      gammaCounter++;
    }
  }

  return fUrgent;
}

void ExN06StackingAction::NewStage()
{
  G4cout << "Number of optical photons produces in this event : "
         << gammaCounter << G4endl;
}

void ExN06StackingAction::PrepareNewEvent()
{ gammaCounter = 0; }
