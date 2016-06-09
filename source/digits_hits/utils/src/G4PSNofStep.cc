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
// $Id: G4PSNofStep.cc,v 1.4 2005/12/13 08:43:51 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// G4PSNofStep
#include "G4PSNofStep.hh"

// (Description)
//   This is a primitive scorer class for scoring number of steps in the
//  Cell.
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
//

G4PSNofStep::G4PSNofStep(G4String name, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1)
{;}

G4PSNofStep::~G4PSNofStep()
{;}

G4bool G4PSNofStep::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4int  index = GetIndex(aStep);
  G4double val = 1.0;
  EvtMap->add(index,val);  
  return TRUE;
}

void G4PSNofStep::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(),GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSNofStep::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSNofStep::clear(){
  EvtMap->clear();
}

void G4PSNofStep::DrawAll()
{;}

void G4PSNofStep::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  num of step: " << *(itr->second)
	   << G4endl;
  }
}

