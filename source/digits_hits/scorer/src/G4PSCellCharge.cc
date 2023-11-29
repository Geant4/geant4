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
// G4PSCellCharge
#include "G4PSCellCharge.hh"
#include "G4Track.hh"

///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring cell charge.
//   The Cell Charge is defined by  a sum of deposited charge inside the cell.
//
// Created: 2006-08-20  Tsukasa ASO
// 2010-07-22   Introduce Unit specification.
//
///////////////////////////////////////////////////////////////////////////////

G4PSCellCharge::G4PSCellCharge(G4String name, G4int depth)
  : G4PSCellCharge(name, "e+", depth) 
{}

G4PSCellCharge::G4PSCellCharge(G4String name, const G4String& unit, G4int depth)
  : G4VPrimitiveScorer(name, depth)
  , HCID(-1)
  , EvtMap(nullptr)
{
  SetUnit(unit);
}

G4bool G4PSCellCharge::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // Enter or First step of primary.
  if(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary ||
     (aStep->GetTrack()->GetParentID() == 0 &&
      aStep->GetTrack()->GetCurrentStepNumber() == 1))
  {
    G4double CellCharge = aStep->GetPreStepPoint()->GetCharge();
    CellCharge *= aStep->GetPreStepPoint()->GetWeight();
    G4int index = GetIndex(aStep);
    EvtMap->add(index, CellCharge);
  }

  // Exit
  if(aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary)
  {
    G4double CellCharge = aStep->GetPreStepPoint()->GetCharge();
    CellCharge *= aStep->GetPreStepPoint()->GetWeight();
    G4int index = GetIndex(aStep);
    CellCharge *= -1.0;
    EvtMap->add(index, CellCharge);
  }

  return true;
}

void G4PSCellCharge::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if(HCID < 0)
    HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, EvtMap);
}

void G4PSCellCharge::clear() { EvtMap->clear(); }

void G4PSCellCharge::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  for(const auto& [copy, charge] : *(EvtMap->GetMap()))
  {
    G4cout << "  copy no.: " << copy
           << "  cell charge : " << *(charge) / GetUnitValue() << " ["
           << GetUnit() << "]" << G4endl;
  }
}

void G4PSCellCharge::SetUnit(const G4String& unit)
{
  CheckAndSetUnit(unit, "Electric charge");
}
