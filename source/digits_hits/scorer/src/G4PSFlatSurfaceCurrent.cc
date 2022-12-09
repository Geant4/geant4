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
// G4PSFlatSurfaceCurrent
#include "G4PSFlatSurfaceCurrent.hh"

#include "G4SystemOfUnits.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
#include "G4GeometryTolerance.hh"
#include "G4VScoreHistFiller.hh"

////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring only Surface Flux.
//  Current version assumes only for G4Box shape.
//
// Surface is defined at the -Z surface.
// Direction                  -Z   +Z
//   0  IN || OUT            ->|<-  |
//   1  IN                   ->|    |
//   2  OUT                    |<-  |
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 17-Nov-2005 T.Aso, Bug fix for area definition.
// 31-Mar-2007 T.Aso, Add option for normalizing by the area.
// 2010-07-22   Introduce Unit specification.
// 2020-10-06   Use G4VPrimitivePlotter and fill 1-D histo of kinetic energy (x)
//              vs. Surface Current * track weight (y)    (Makoto Asai)
//
///////////////////////////////////////////////////////////////////////////////

G4PSFlatSurfaceCurrent::G4PSFlatSurfaceCurrent(G4String name, G4int direction,
                                               G4int depth)
  : G4PSFlatSurfaceCurrent(name, direction, "percm2", depth) 
{}

G4PSFlatSurfaceCurrent::G4PSFlatSurfaceCurrent(G4String name, G4int direction,
                                               const G4String& unit,
                                               G4int depth)
  : G4VPrimitivePlotter(name, depth)
  , HCID(-1)
  , fDirection(direction)
  , EvtMap(nullptr)
  , weighted(true)
  , divideByArea(true)
{
  DefineUnitAndCategory();
  SetUnit(unit);
}

G4bool G4PSFlatSurfaceCurrent::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4StepPoint* preStep             = aStep->GetPreStepPoint();
  G4VPhysicalVolume* physVol       = preStep->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid* solid                  = nullptr;
  if(physParam != nullptr)
  {  // for parameterized volume
    G4int idx =
      ((G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable()))
        ->GetReplicaNumber(indexDepth);
    solid = physParam->ComputeSolid(idx, physVol);
    solid->ComputeDimensions(physParam, idx, physVol);
  }
  else
  {  // for ordinary volume
    solid = physVol->GetLogicalVolume()->GetSolid();
  }

  auto  boxSolid = (G4Box*) (solid);

  G4int dirFlag = IsSelectedSurface(aStep, boxSolid);
  if(dirFlag > 0)
  {
    if(fDirection == fCurrent_InOut || fDirection == dirFlag)
    {
      G4int index                    = GetIndex(aStep);
      G4TouchableHandle theTouchable = preStep->GetTouchableHandle();
      G4double current               = 1.0;
      if(weighted)
        current = preStep->GetWeight();  // Current (Particle Weight)
      if(divideByArea)
      {
        G4double square =
          4. * boxSolid->GetXHalfLength() * boxSolid->GetYHalfLength();
        current = current / square;  // Normalized by Area
      }
      EvtMap->add(index, current);

      if(!hitIDMap.empty() && hitIDMap.find(index) != hitIDMap.cend())
      {
        auto filler = G4VScoreHistFiller::Instance();
        if(filler == nullptr)
        {
          G4Exception("G4PSFlatSurfaceCurrent::ProcessHits", "SCORER0123",
                      JustWarning,
                      "G4TScoreHistFiller is not instantiated!! Histogram is "
                      "not filled.");
        }
        else
        {
          filler->FillH1(hitIDMap[index], preStep->GetKineticEnergy(), current);
        }
      }
    }
  }

  return true;
}

G4int G4PSFlatSurfaceCurrent::IsSelectedSurface(G4Step* aStep, G4Box* boxSolid)
{
  G4TouchableHandle theTouchable =
    aStep->GetPreStepPoint()->GetTouchableHandle();
  G4double kCarTolerance =
    G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  if(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary)
  {
    // Entering Geometry
    G4ThreeVector stppos1 = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector localpos1 =
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos1);
    if(std::fabs(localpos1.z() + boxSolid->GetZHalfLength()) < kCarTolerance)
    {
      return fCurrent_In;
    }
  }

  if(aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary)
  {
    // Exiting Geometry
    G4ThreeVector stppos2 = aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector localpos2 =
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos2);
    if(std::fabs(localpos2.z() + boxSolid->GetZHalfLength()) < kCarTolerance)
    {
      return fCurrent_Out;
    }
  }

  return -1;
}

void G4PSFlatSurfaceCurrent::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if(HCID < 0)
    HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*) EvtMap);
}

void G4PSFlatSurfaceCurrent::clear() { EvtMap->clear(); }

void G4PSFlatSurfaceCurrent::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  for(const auto& [copy, current] : *(EvtMap->GetMap()))
  {
    G4cout << "  copy no.: " << copy << " current : ";
    if(divideByArea)
    {
      G4cout << *(current) / GetUnitValue() << " [" << GetUnit() << "]";
    }
    else
    {
      G4cout << *(current) / GetUnitValue() << " [tracks]";
    }
    G4cout << G4endl;
  }
}

void G4PSFlatSurfaceCurrent::SetUnit(const G4String& unit)
{
  if(divideByArea)
  {
    CheckAndSetUnit(unit, "Per Unit Surface");
  }
  else
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
      G4Exception("G4PSFlatSurfaceCurrent::SetUnit", "DetPS0007", JustWarning,
                  msg);
    }
  }
}

void G4PSFlatSurfaceCurrent::DefineUnitAndCategory()
{
  // Per Unit Surface
  new G4UnitDefinition("percentimeter2", "percm2", "Per Unit Surface",
                       (1. / cm2));
  new G4UnitDefinition("permillimeter2", "permm2", "Per Unit Surface",
                       (1. / mm2));
  new G4UnitDefinition("permeter2", "perm2", "Per Unit Surface", (1. / m2));
}
