//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4AdjointStackingAction class implementation
//
// Author: L. Desorgher, SpaceIT GmbH - April 2008
// Contract: ESA contract 21435/08/NL/AT
// Customer: ESA/ESTEC
// --------------------------------------------------------------------

#include "G4AdjointStackingAction.hh"
#include "G4AdjointTrackingAction.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4StackManager.hh"

G4AdjointStackingAction::
G4AdjointStackingAction(G4AdjointTrackingAction* anAction)
{
  theAdjointTrackingAction = anAction;
}

// --------------------------------------------------------------------
//
G4ClassificationOfNewTrack
G4AdjointStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fUrgent;
  G4String partType = aTrack->GetParticleDefinition()->GetParticleType();
  adjoint_mode = G4StrUtil::contains(partType, "adjoint");
  if (!adjoint_mode )
  {
    if (!reclassification_stage)
    {
      classification = fWaiting;
    }
    else  // need to check if forwrad tracking can be continued use of
    {
      if (theAdjointTrackingAction->GetNbOfAdointTracksReachingTheExternalSurface()>0)
      {
        if (theFwdStackingAction != nullptr)
        {
          classification =  theFwdStackingAction->ClassifyNewTrack(aTrack);
        }
      }
      else
      {
        classification = fKill;
      }
    }
  }
  else if (theUserAdjointStackingAction != nullptr)
  {
    classification = theUserAdjointStackingAction->ClassifyNewTrack(aTrack);
  }
  return classification;
}

// --------------------------------------------------------------------
//
void G4AdjointStackingAction::NewStage()
{
  reclassification_stage = true;
  if (first_reclassification_stage)
  {
    if (theUserAdjointStackingAction != nullptr)
    {
      theUserAdjointStackingAction->NewStage();
    }
    stackManager->ReClassify();
  }
  else if (theFwdStackingAction != nullptr) theFwdStackingAction->NewStage();
  {
    first_reclassification_stage = false;
  }
}

// --------------------------------------------------------------------
//    
void G4AdjointStackingAction::PrepareNewEvent()
{
  reclassification_stage = false;
  first_reclassification_stage = true;
  if (theUserAdjointStackingAction != nullptr)
  {
    theUserAdjointStackingAction->PrepareNewEvent();
  }
}
