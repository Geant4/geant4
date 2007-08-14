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
// $Id: G4PSTermination.cc,v 1.1 2007-08-14 21:23:52 taso Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4PSTermination
#include "G4PSTermination.hh"

///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring only terminated traks inside
//  the cell.
//
// Created: 2007-02-02  Tsukasa ASO, Akinori Kimura.
//
///////////////////////////////////////////////////////////////////////////////

G4PSTermination::G4PSTermination(G4String name, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1),weighted(false)
{;}

G4PSTermination::~G4PSTermination()
{;}

G4bool G4PSTermination::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  if ( aStep->GetTrack()->GetTrackStatus() != fStopAndKill ) return FALSE;

  G4int  index = GetIndex(aStep);
  G4double val = 1.0;
  if(weighted) val *= aStep->GetPreStepPoint()->GetWeight();
  EvtMap->add(index,val);  
  return TRUE;
}

void G4PSTermination::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(),GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSTermination::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSTermination::clear(){
  EvtMap->clear();
}

void G4PSTermination::DrawAll()
{;}

void G4PSTermination::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  collisions: " << *(itr->second)
	   << G4endl;
  }
}

