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
// G4PSPassageCellCurrent
#include "G4PSPassageCellCurrent.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4UnitsTable.hh"
#include "G4VScoreHistFiller.hh"

////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring number of tracks,
//   where only tracks passing through the geometry are taken
//  into account.
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 2010-07-22   Introduce Unit specification.
// 2010-07-22   Add weighted option
// 2020-10-06   Use G4VPrimitivePlotter and fill 1-D histo of kinetic energy (x)
//              vs. current * track weight (y)               (Makoto Asai)
//
///////////////////////////////////////////////////////////////////////////////

G4PSPassageCellCurrent::G4PSPassageCellCurrent(G4String name, G4int depth)
  : G4VPrimitivePlotter(name, depth)
{
  SetUnit("");
}

G4bool G4PSPassageCellCurrent::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  if(IsPassed(aStep))
  {
    fCurrent = 1.;
    if(weighted)
      fCurrent = aStep->GetPreStepPoint()->GetWeight();
    G4int index = GetIndex(aStep);
    EvtMap->add(index, fCurrent);

    if(!hitIDMap.empty() && hitIDMap.find(index) != hitIDMap.cend())
    {
      auto filler = G4VScoreHistFiller::Instance();
      if(filler == nullptr)
      {
        G4Exception(
          "G4PSVolumeFlux::ProcessHits", "SCORER0123", JustWarning,
          "G4TScoreHistFiller is not instantiated!! Histogram is not filled.");
      }
      else
      {
        filler->FillH1(hitIDMap[index],
                       aStep->GetPreStepPoint()->GetKineticEnergy(), fCurrent);
      }
    }
  }

  return true;
}

G4bool G4PSPassageCellCurrent::IsPassed(G4Step* aStep)
{
  G4bool Passed = false;

  G4bool IsEnter = aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary;
  G4bool IsExit  = aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary;

  G4int trkid = aStep->GetTrack()->GetTrackID();

  if(IsEnter && IsExit)
  {  // Passed at one step
    Passed = true;
  }
  else if(IsEnter)
  {                         // Enter a new geometry
    fCurrentTrkID = trkid;  // Resetting the current track.
  }
  else if(IsExit)
  {  // Exit a current geometry
    if(fCurrentTrkID == trkid)
    {
      Passed = true;  // if the track is same as entered.
    }
  }
  else
  {  // Inside geometry
    if(fCurrentTrkID == trkid)
    {  // Adding the track length to current one ,
    }
  }
  return Passed;
}

void G4PSPassageCellCurrent::Initialize(G4HCofThisEvent* HCE)
{
  fCurrentTrkID = -1;

  EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if(HCID < 0)
    HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, EvtMap);
}

void G4PSPassageCellCurrent::clear() { EvtMap->clear(); }

void G4PSPassageCellCurrent::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  for(const auto& [copy, current] : *(EvtMap->GetMap()))
  {
    G4cout << "  copy no.: " << copy
           << "  cell current : " << *(current) << " [tracks] " << G4endl;
  }
}

void G4PSPassageCellCurrent::SetUnit(const G4String& unit)
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
    G4Exception("G4PSPassageCellCurrent::SetUnit", "DetPS0012", JustWarning,
                msg);
  }
}
