// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02StackingAction.cc,v 1.2 1999/11/29 18:33:29 morita Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

#include "PersEx02StackingAction.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"

PersEx02StackingAction::PersEx02StackingAction()
{;}

PersEx02StackingAction::~PersEx02StackingAction()
{;}

G4ClassificationOfNewTrack 
PersEx02StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fKill;
  if(aTrack->GetParentID()==0)
  { classification = fUrgent; }
  return classification;
}

void PersEx02StackingAction::NewStage()
{;}
    
void PersEx02StackingAction::PrepareNewEvent()
{;}


