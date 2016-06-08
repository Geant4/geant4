
#include "ExE02StackingAction.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"

ExE02StackingAction::ExE02StackingAction()
{;}

ExE02StackingAction::~ExE02StackingAction()
{;}

G4ClassificationOfNewTrack 
ExE02StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fKill;
  if(aTrack->GetParentID()==0)
  { classification = fUrgent; }
  return classification;
}

void ExE02StackingAction::NewStage()
{;}
    
void ExE02StackingAction::PrepareNewEvent()
{;}


