// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserStackingAction.cc,v 1.2.8.1 1999/12/07 20:47:55 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

#include "G4UserStackingAction.hh"
#include "G4Track.hh"
#include "G4ios.hh"

G4UserStackingAction::G4UserStackingAction()
{;}

G4UserStackingAction::~G4UserStackingAction()
{;}

G4ClassificationOfNewTrack G4UserStackingAction::ClassifyNewTrack
(const G4Track* aTrack)
{
  return fUrgent;
}

void G4UserStackingAction::NewStage()
{;}

void G4UserStackingAction::PrepareNewEvent()
{;}


