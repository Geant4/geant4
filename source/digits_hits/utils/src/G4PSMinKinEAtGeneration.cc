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
// $Id: G4PSMinKinEAtGeneration.cc,v 1.4 2005/12/13 08:43:45 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// G4PSMinKinEAtGeneration
#include "G4PSMinKinEAtGeneration.hh"
#include "G4UnitsTable.hh"

// (Description)
//   This is a primitive scorer class for scoring minimum energy of
//  newly produced particles in the geometry.
//
// Created: 2005-11-17  Tsukasa ASO, Akinori Kimura.
//

G4PSMinKinEAtGeneration::G4PSMinKinEAtGeneration(G4String name, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1)
{;}

G4PSMinKinEAtGeneration::~G4PSMinKinEAtGeneration()
{;}

G4bool G4PSMinKinEAtGeneration::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  // ===============
  //  First Step, 
  //  Confirm this is a newly produced secondary.
  //  
  //- check for newly produced particle. e.g. Step number is 1.
  if ( aStep->GetTrack()->GetCurrentStepNumber() != 1) return FALSE;
  //- check for this is not a primary particle. e.g. ParentID != 0 .
  if ( aStep->GetTrack()->GetParentID() == 0 ) return FALSE;

  //===============================================
  //- This is a newly produced secondary particle.
  //===============================================

  // ===============
  //  Second Step,
  //  Confirm this track has lower energy than previous one.
  //

  // -Kinetic energy of this particle at the starting point.
  G4double kinetic = aStep->GetPreStepPoint()->GetKineticEnergy();

  // -Stored value in the current HitsMap.
  G4int  index = GetIndex(aStep);
  G4double*   mapValue=  ((*EvtMap)[index]);

  // -
  // If mapValue exits (e.g not NULL ), compare it with 
  // current track's kinetic energy.
  if ( mapValue && ( kinetic > *mapValue ) ) return FALSE;

  // -
  // Current Track is a newly produced secondary and has lower
  // kinetic energy than previous one in this geometry.
  //
  EvtMap->set(index, kinetic);
  return TRUE;
}

void G4PSMinKinEAtGeneration::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(),GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSMinKinEAtGeneration::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSMinKinEAtGeneration::clear(){
  EvtMap->clear();
}

void G4PSMinKinEAtGeneration::DrawAll()
{;}

void G4PSMinKinEAtGeneration::PrintAll()
{
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  num of step: " << G4BestUnit(*(itr->second),"Energy")
	   << G4endl;
  }
}

