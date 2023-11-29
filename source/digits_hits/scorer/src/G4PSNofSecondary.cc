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
// G4PSNofSecondary
#include "G4PSNofSecondary.hh"
#include "G4VScoreHistFiller.hh"

// (Description)
//   This is a primitive scorer class for scoring number of secondaries
// in the Cell.
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// Modify:  2011-09-09  T.Aso modify comment in PrintAll().
//          2020-10-06   Use G4VPrimitivePlotter and fill 1-D histo of kinetic
//                       energy of the secondary in MeV (x) vs. track weight (y)
//                       (Makoto Asai)
//

G4PSNofSecondary::G4PSNofSecondary(G4String name, G4int depth)
  : G4VPrimitivePlotter(name, depth)
{}

G4bool G4PSNofSecondary::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  //- check for newly produced particle. e.g. first step.
  if(aStep->GetTrack()->GetCurrentStepNumber() != 1)
    return false;
  //- check for this is not a primary particle. e.g. ParentID > 0 .
  if(aStep->GetTrack()->GetParentID() == 0)
    return false;
  //- check the particle if the partifle definition is given.
  if((particleDef != nullptr) && particleDef != aStep->GetTrack()->GetDefinition())
    return false;
  //
  //- This is a newly produced secondary particle.
  G4int index     = GetIndex(aStep);
  G4double weight = 1.0;
  if(weighted)
    weight *= aStep->GetPreStepPoint()->GetWeight();
  EvtMap->add(index, weight);

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
                     aStep->GetPreStepPoint()->GetKineticEnergy(), weight);
    }
  }

  return true;
}

void G4PSNofSecondary::SetParticle(const G4String& particleName)
{
  G4ParticleDefinition* pd =
    G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  if(pd == nullptr)
  {
    G4String msg = "Particle <";
    msg += particleName;
    msg += "> not found.";
    G4Exception("G4PSNofSecondary::SetParticle", "DetPS0101", FatalException,
                msg);
  }
  particleDef = pd;
}

void G4PSNofSecondary::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if(HCID < 0)
  {
    HCID = GetCollectionID(0);
  }
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*) EvtMap);
}

void G4PSNofSecondary::clear() { EvtMap->clear(); }

void G4PSNofSecondary::PrintAll()
{
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  for(const auto& [copy, secondaries] : *(EvtMap->GetMap()))
  {
    G4cout << "  copy no.: " << copy
           << "  num of secondaries: " << *(secondaries) / GetUnitValue()
           << G4endl;
  }
}

void G4PSNofSecondary::SetUnit(const G4String& unit)
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
    G4Exception("G4PSNofSecondary::SetUnit", "DetPS0010", JustWarning, msg);
  }
}
