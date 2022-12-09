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
// G4PSTrackLength
#include "G4PSTrackLength.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring sum of track length.
//
//
// Created: 2007-02-02  Tsukasa ASO, Akinori Kimura.
//          2010-07-22  Introduce Unit specification.
//          2011-09-09  Modify comment in PrintAll().
//
///////////////////////////////////////////////////////////////////////////////

G4PSTrackLength::G4PSTrackLength(G4String name, G4int depth)
  : G4PSTrackLength(name, "mm", depth)
{}

G4PSTrackLength::G4PSTrackLength(G4String name, const G4String& unit,
                                 G4int depth)
  : G4VPrimitiveScorer(name, depth)
  , HCID(-1)
  , EvtMap(nullptr)
  , weighted(false)
  , multiplyKinE(false)
  , divideByVelocity(false)
{
  DefineUnitAndCategory();
  SetUnit(unit);
}

void G4PSTrackLength::MultiplyKineticEnergy(G4bool flg)
{
  multiplyKinE = flg;
  // Default unit is set according to flags.
  SetUnit("");
}

void G4PSTrackLength::DivideByVelocity(G4bool flg)
{
  divideByVelocity = flg;
  // Default unit is set according to flags.
  SetUnit("");
}

G4bool G4PSTrackLength::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double trklength = aStep->GetStepLength();
  if(trklength == 0.)
    return false;
  if(weighted)
    trklength *= aStep->GetPreStepPoint()->GetWeight();
  if(multiplyKinE)
    trklength *= aStep->GetPreStepPoint()->GetKineticEnergy();
  if(divideByVelocity)
    trklength /= aStep->GetPreStepPoint()->GetVelocity();
  G4int index = GetIndex(aStep);
  EvtMap->add(index, trklength);
  return true;
}

void G4PSTrackLength::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if(HCID < 0)
  {
    HCID = GetCollectionID(0);
  }
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*) EvtMap);
}

void G4PSTrackLength::clear() { EvtMap->clear(); }

void G4PSTrackLength::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  for(const auto& [copy, length] : *(EvtMap->GetMap()))
  {
    G4cout << "  copy no.: " << copy;
    if(multiplyKinE)
    {
      if(divideByVelocity)
        G4cout << " EnergyFlux: ";
      else
        G4cout << " EnergyFlow: ";
    }
    else
    {
      if(divideByVelocity)
        G4cout << " Time: ";
      else
        G4cout << " Length: ";
    }
    G4cout << *(length) / GetUnitValue() << " [" << GetUnit() << "]";
    G4cout << G4endl;
  }
}

void G4PSTrackLength::SetUnit(const G4String& unit)
{
  if(multiplyKinE)
  {
    if(divideByVelocity)
    {
      if(unit.empty())
      {
        CheckAndSetUnit("MeV_second", "EnergyFlux");
      }
      else
      {
        CheckAndSetUnit(unit, "EnergyFlux");
      }
    }
    else
    {
      if(unit.empty())
      {
        CheckAndSetUnit("MeV_mm", "EnergyFlow");
      }
      else
      {
        CheckAndSetUnit(unit, "EnergyFlow");
      }
    }
  }
  else
  {
    if(divideByVelocity)
    {
      if(unit.empty())
      {
        CheckAndSetUnit("second", "Time");
      }
      else
      {
        CheckAndSetUnit(unit, "Time");
      }
    }
    else
    {
      if(unit.empty())
      {
        CheckAndSetUnit("mm", "Length");
      }
      else
      {
        CheckAndSetUnit(unit, "Length");
      }
    }
  }
}

void G4PSTrackLength::DefineUnitAndCategory()
{
  // EnergyFlux
  new G4UnitDefinition("eV_second", "eV_s", "EnergyFlux", (eV * second));
  new G4UnitDefinition("keV_second", "keV_s", "EnergyFlux", (keV * second));
  new G4UnitDefinition("MeV_second", "MeV_s", "EnergyFlux", (MeV * second));
  new G4UnitDefinition("eV_millisecond", "eV_ms", "EnergyFlux", (eV * ms));
  new G4UnitDefinition("keV_millisecond", "keV_ms", "EnergyFlux", (keV * ms));
  new G4UnitDefinition("MeV_millisecond", "MeV_ms", "EnergyFlux", (MeV * ms));
  // EnergyFlow
  new G4UnitDefinition("eV_millimeter", "eV_mm", "EnergyFlow", (eV * mm));
  new G4UnitDefinition("keV_millimeter", "keV_mm", "EnergyFlow", (keV * mm));
  new G4UnitDefinition("MeV_millimeter", "MeV_mm", "EnergyFlow", (MeV * mm));
  new G4UnitDefinition("eV_centimeter", "eV_cm", "EnergyFlow", (eV * cm));
  new G4UnitDefinition("keV_centimeter", "keV_cm", "EnergyFlow", (keV * cm));
  new G4UnitDefinition("MeV_centimeter", "MeV_cm", "EnergyFlow", (MeV * cm));
  new G4UnitDefinition("eV_meter", "eV_m", "EnergyFlow", (eV * m));
  new G4UnitDefinition("keV_meter", "keV_m", "EnergyFlow", (keV * m));
  new G4UnitDefinition("MeV_meter", "MeV_m", "EnergyFlow", (MeV * m));
}
