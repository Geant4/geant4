
#include "PersEx01StackingAction.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"

PersEx01StackingAction::PersEx01StackingAction()
{;}

PersEx01StackingAction::~PersEx01StackingAction()
{;}

G4ClassificationOfNewTrack 
PersEx01StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fKill;
  if(aTrack->GetParentID()==0)
  { classification = fUrgent; }
  return classification;
}

void PersEx01StackingAction::NewStage()
{;}
    
void PersEx01StackingAction::PrepareNewEvent()
{;}


