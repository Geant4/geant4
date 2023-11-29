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
// Class G4PSVolumeFlux
//
// Class description:
//  Scorer that scores number of tracks coming into the associated volume.
//  Optionally number of tracks can be divided by the surface area to
//  score the volume current and divided by cos(theta) for volume flux
//  where theta is the incident angle
//
//  - Created   M. Asai, Sept. 2020
//
//
#include "G4PSVolumeFlux.hh"

#include "G4SystemOfUnits.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
#include "G4GeometryTolerance.hh"
#include "G4VScoreHistFiller.hh"

G4PSVolumeFlux::G4PSVolumeFlux(G4String name, G4int direction, G4int depth)
  : G4VPrimitivePlotter(name, depth)
  , fDirection(direction)
{}

G4bool G4PSVolumeFlux::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4StepPoint* preStepPoint  = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  G4StepPoint* thisStepPoint = nullptr;
  if(fDirection == 1)  // Score only the inward particle
  {
    if(preStepPoint->GetStepStatus() == fGeomBoundary)
    {
      thisStepPoint = preStepPoint;
    }
    else
    {
      return false;
    }
  }
  else if(fDirection == 2)  // Score only the outward particle
  {
    if(postStepPoint->GetStepStatus() == fGeomBoundary)
    {
      thisStepPoint = postStepPoint;
    }
    else
    {
      return false;
    }
  }

  G4double flux = preStepPoint->GetWeight();

  if(divare || divcos)
  {
    G4VPhysicalVolume* physVol       = preStepPoint->GetPhysicalVolume();
    G4VPVParameterisation* physParam = physVol->GetParameterisation();
    G4VSolid* solid                  = nullptr;
    if(physParam != nullptr)
    {  // for parameterized volume
      auto idx = ((G4TouchableHistory*) (preStepPoint->GetTouchable()))
                   ->GetReplicaNumber(indexDepth);
      solid = physParam->ComputeSolid(idx, physVol);
      solid->ComputeDimensions(physParam, idx, physVol);
    }
    else
    {  // for ordinary volume
      solid = physVol->GetLogicalVolume()->GetSolid();
    }

    if(divare)
    {
      flux /= solid->GetSurfaceArea();
    }

    if(divcos)
    {
      G4TouchableHandle theTouchable = thisStepPoint->GetTouchableHandle();
      G4ThreeVector pdirection       = thisStepPoint->GetMomentumDirection();
      G4ThreeVector localdir =
        theTouchable->GetHistory()->GetTopTransform().TransformAxis(pdirection);
      G4ThreeVector globalPos = thisStepPoint->GetPosition();
      G4ThreeVector localPos =
        theTouchable->GetHistory()->GetTopTransform().TransformPoint(globalPos);
      G4ThreeVector surfNormal = solid->SurfaceNormal(localPos);
      G4double cosT            = surfNormal.cosTheta(localdir);
      if(cosT != 0.)
        flux /= std::abs(cosT);
    }
  }

  G4int index = GetIndex(aStep);
  EvtMap->add(index, flux);

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
      filler->FillH1(hitIDMap[index], thisStepPoint->GetKineticEnergy(), flux);
    }
  }

  return true;
}

void G4PSVolumeFlux::Initialize(G4HCofThisEvent* HCE)
{
  if(HCID < 0)
    HCID = GetCollectionID(0);
  EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
                                    GetName());
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*) EvtMap);
}

void G4PSVolumeFlux::clear() { EvtMap->clear(); }

void G4PSVolumeFlux::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer" << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  for(const auto& [copy, flux] : *(EvtMap->GetMap()))
  {
    G4cout << "  copy no.: " << copy << "  flux  : " << *(flux)
           << G4endl;
  }
}
