// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserStackingAction.cc,v 1.3 1999-12-15 14:49:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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


