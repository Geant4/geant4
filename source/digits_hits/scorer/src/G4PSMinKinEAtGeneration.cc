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
//
// $Id: G4PSMinKinEAtGeneration.cc 81087 2014-05-20 15:44:27Z gcosmo $
//
// G4PSMinKinEAtGeneration
#include "G4PSMinKinEAtGeneration.hh"
#include "G4UnitsTable.hh"

// (Description)
//   This is a primitive scorer class for scoring minimum energy of
//  newly produced particles in the geometry.
//
// Created: 2005-11-17  Tsukasa ASO, Akinori Kimura.
// 2010-07-22   Introduce Unit specification.
// 2011-09-09   Modify comment in PrintAll().  T.Aso
//

G4PSMinKinEAtGeneration::G4PSMinKinEAtGeneration(G4String name, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0)
{
    SetUnit("MeV");
}

G4PSMinKinEAtGeneration::G4PSMinKinEAtGeneration(G4String name, 
						 const G4String& unit,
						 G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0)
{
    SetUnit(unit);
}

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
	   << "  energy: " << *(itr->second)/GetUnitValue()
	   << " ["<<GetUnit()<<"]"
	   << G4endl;
  }
}

void G4PSMinKinEAtGeneration::SetUnit(const G4String& unit)
{
    CheckAndSetUnit(unit,"Energy");
}


