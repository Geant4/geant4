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
// $Id: G4PSNofSecondary.cc,v 1.2 2005/12/13 08:43:49 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// G4PSNofSecondary
#include "G4PSNofSecondary.hh"

// (Description)
//   This is a primitive scorer class for scoring number of steps in the
//  Cell.
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
//

G4PSNofSecondary::G4PSNofSecondary(G4String name, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1)
{;}

G4PSNofSecondary::~G4PSNofSecondary()
{;}

G4bool G4PSNofSecondary::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  //- check for newly produced particle. e.g. first step.
  if ( aStep->GetTrack()->GetCurrentStepNumber() != 1) return FALSE;
  //- check for this is not a primary particle. e.g. ParentID > 0 .
  if ( aStep->GetTrack()->GetParentID() == 0 ) return FALSE;
  //
  //- This is a newly produced secondary particle.
  G4int  index = GetIndex(aStep);
  G4double weight = aStep->GetPreStepPoint()->GetWeight();
  EvtMap->add(index,weight);  
  return TRUE;
}

void G4PSNofSecondary::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(),GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSNofSecondary::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSNofSecondary::clear(){
  EvtMap->clear();
}

void G4PSNofSecondary::DrawAll()
{;}

void G4PSNofSecondary::PrintAll()
{
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  num of step: " << *(itr->second)
	   << G4endl;
  }
}

