//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: PersEx02StackingAction.cc,v 1.3 2001/07/11 09:58:15 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
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


