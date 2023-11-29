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
// G4PSTrackCounter
#include "G4PSTrackCounter.hh"
#include "G4UnitsTable.hh"
#include "G4VScoreHistFiller.hh"

///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring number of tracks in a cell.
//
// Created: 2007-02-02  Tsukasa ASO, Akinori Kimura.
// 2010-07-22   Introduce Unit specification.
//
///////////////////////////////////////////////////////////////////////////////

G4PSTrackCounter::G4PSTrackCounter(G4String name, G4int direction, G4int depth)
  : G4VPrimitivePlotter(name, depth)
  , fDirection(direction)
{
  SetUnit("");
}

G4bool G4PSTrackCounter::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4StepPoint* posStep = aStep->GetPostStepPoint();
  G4bool IsEnter       = preStep->GetStepStatus() == fGeomBoundary;
  G4bool IsExit        = posStep->GetStepStatus() == fGeomBoundary;

  // Regular : count in prestep volume.
  G4int index = GetIndex(aStep);

  //  G4cout << " trk " << aStep->GetTrack()->GetTrackID()
  //	 << " index " << index << " In " << IsEnter << " Out " <<IsExit
  //	 << " PreStatus " << preStep->GetStepStatus()
  //	 << " PostStatus " << posStep->GetStepStatus()
  //	 << " w " << aStep->GetPreStepPoint()->GetWeight()
  // <<    G4endl;

  G4bool flag = false;

  if(IsEnter && fDirection == fCurrent_In)
    flag = true;
  else if(IsExit && fDirection == fCurrent_Out)
    flag = true;
  else if((IsExit || IsEnter) && fDirection == fCurrent_InOut)
    flag = true;

  if(flag)
  {
    // G4cout << " Count " << G4endl;
    G4double val = 1.0;
    if(weighted)
      val *= aStep->GetPreStepPoint()->GetWeight();
    EvtMap->add(index, val);

    if(!hitIDMap.empty() && hitIDMap.find(index) != hitIDMap.cend())
    {
      auto filler = G4VScoreHistFiller::Instance();
      if(filler == nullptr)
      {
        G4Exception(
          "G4PSTrackCounter::ProcessHits", "SCORER0123", JustWarning,
          "G4TScoreHistFiller is not instantiated!! Histogram is not filled.");
      }
      else
      {
        filler->FillH1(hitIDMap[index],
                       aStep->GetPreStepPoint()->GetKineticEnergy(), val);
      }
    }
  }

  return true;
}

void G4PSTrackCounter::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if(HCID < 0)
  {
    HCID = GetCollectionID(0);
  }
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*) EvtMap);
  // G4cout <<   "----- PS Initialize " << G4endl;
}

void G4PSTrackCounter::clear() { EvtMap->clear(); }

void G4PSTrackCounter::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  for(const auto& [copy, count] : *(EvtMap->GetMap()))
  {
    G4cout << "  copy no.: " << copy
           << "  track count: " << *(count) << " [tracks] " << G4endl;
  }
}

void G4PSTrackCounter::SetUnit(const G4String& unit)
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
    G4Exception("G4PSTrackCounter::SetUnit", "DetPS0018", JustWarning, msg);
  }
}
