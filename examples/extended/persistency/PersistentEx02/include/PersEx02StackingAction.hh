// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02StackingAction.hh,v 1.2 1999/11/29 18:23:33 morita Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#ifndef PersEx02StackingAction_H
#define PersEx02StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"

class G4Track;

class PersEx02StackingAction : public G4UserStackingAction
{
  public:
    PersEx02StackingAction();
    virtual ~PersEx02StackingAction();

  public:
   virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
   virtual void NewStage();
   virtual void PrepareNewEvent();

};

#endif

