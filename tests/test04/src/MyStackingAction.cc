
#include "MyStackingAction.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"

MyStackingAction::MyStackingAction()
{;}

MyStackingAction::~MyStackingAction()
{;}

G4ClassificationOfNewTrack 
 MyStackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  G4ClassificationOfNewTrack classification = fKill;
  if(aTrack->GetParentID()==0)
  { classification = fUrgent; }
  return classification;
}

void MyStackingAction::NewStage()
{;}
    
void MyStackingAction::PrepareNewEvent()
{;}


