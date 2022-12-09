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
//
// G4PSPopulation
#include "G4PSPopulation.hh"

///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring population in a cell.
//
// Created: 2007-02-02  Tsukasa ASO, Akinori Kimura.
// 2010-07-22   Introduce Unit specification.
//
///////////////////////////////////////////////////////////////////////////////

G4PSPopulation::G4PSPopulation(G4String name, G4int depth)
  : G4VPrimitiveScorer(name, depth)
{
  SetUnit("");
}

G4bool G4PSPopulation::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4int index         = GetIndex(aStep);
  G4TrackLogger& tlog = fCellTrackLogger[index];
  if(tlog.FirstEnterance(aStep->GetTrack()->GetTrackID()))
  {
    G4double val = 1.0;
    if(weighted)
      val *= aStep->GetPreStepPoint()->GetWeight();
    EvtMap->add(index, val);
  }

  return true;
}

void G4PSPopulation::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if(HCID < 0)
  {
    HCID = GetCollectionID(0);
  }
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*) EvtMap);
}

void G4PSPopulation::EndOfEvent(G4HCofThisEvent*) { fCellTrackLogger.clear(); }

void G4PSPopulation::clear()
{
  EvtMap->clear();
  fCellTrackLogger.clear();
}

void G4PSPopulation::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  for(const auto& [copy, population] : *(EvtMap->GetMap()))
  {
    G4cout << "  copy no.: " << copy
           << "  population: " << *(population) / GetUnitValue() << " [tracks]"
           << G4endl;
  }
}

void G4PSPopulation::SetUnit(const G4String& unit)
{
  if(unit.empty())
  {
    unitName  = unit;
    unitValue = 1.0;
  }
  else
  {
    G4String msg = "Invalid unit [" + unit + "] (Current  unit is [" +
                   GetUnit() + "] ) for " + GetName();
    G4Exception("G4PSPopulation::SetUnit", "DetPS0014", JustWarning, msg);
  }
}
