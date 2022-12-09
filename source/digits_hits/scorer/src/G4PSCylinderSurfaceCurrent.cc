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
// G4PSCylinderSurfaceCurrent
#include "G4PSCylinderSurfaceCurrent.hh"

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
//   This is a primitive scorer class for scoring only Surface Current.
//  Current version assumes only for G4Tubs shape.
//
// Surface is defined at the inner surface of the tube.
// Direction                   R    R+dR
//   0  IN || OUT            ->|<-  |
//   1  IN                   ->|    |
//   2  OUT                    |<-  |
//
// Created: 2007-03-21  Tsukasa ASO
// 2010-07-22   Introduce Unit specification.
// 2020-10-06   Use G4VPrimitivePlotter and fill 1-D histo of kinetic energy (x)
//              vs. surface current * track weight (y)             (Makoto Asai)
//
///////////////////////////////////////////////////////////////////////////////

G4PSCylinderSurfaceCurrent::G4PSCylinderSurfaceCurrent(G4String name,
                                                       G4int direction,
                                                       G4int depth)
  : G4PSCylinderSurfaceCurrent(name, direction, "percm2", depth) 
{}

G4PSCylinderSurfaceCurrent::G4PSCylinderSurfaceCurrent(G4String name,
                                                       G4int direction,
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

G4bool G4PSCylinderSurfaceCurrent::ProcessHits(G4Step* aStep,
                                               G4TouchableHistory*)
{
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4VSolid* solid      = ComputeCurrentSolid(aStep);
  auto  tubsSolid    = static_cast<G4Tubs*>(solid);

  G4int dirFlag = IsSelectedSurface(aStep, tubsSolid);
  // G4cout << " pos " << preStep->GetPosition() <<" dirFlag " << G4endl;
  if(dirFlag > 0)
  {
    if(fDirection == fCurrent_InOut || fDirection == dirFlag)
    {
      G4TouchableHandle theTouchable = preStep->GetTouchableHandle();
      //
      G4double current = 1.0;
      if(weighted)
        current = preStep->GetWeight();  // Current (Particle Weight)
      //
      if(divideByArea)
      {
        G4double square = 2. * tubsSolid->GetZHalfLength() *
                          tubsSolid->GetInnerRadius() *
                          tubsSolid->GetDeltaPhiAngle() / radian;
        current = current / square;  // Current normalized by Area
      }

      G4int index = GetIndex(aStep);
      EvtMap->add(index, current);

      if(!hitIDMap.empty() && hitIDMap.find(index) != hitIDMap.cend())
      {
        auto filler = G4VScoreHistFiller::Instance();
        if(filler == nullptr)
        {
          G4Exception("G4PSCylinderSurfaceCurrent::ProcessHits", "SCORER0123",
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

G4int G4PSCylinderSurfaceCurrent::IsSelectedSurface(G4Step* aStep,
                                                    G4Tubs* tubsSolid)
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
    if(std::fabs(localpos1.z()) > tubsSolid->GetZHalfLength())
      return -1;
    G4double localR2 =
      localpos1.x() * localpos1.x() + localpos1.y() * localpos1.y();
    G4double InsideRadius = tubsSolid->GetInnerRadius();
    if(localR2 >
         (InsideRadius - kCarTolerance) * (InsideRadius - kCarTolerance) &&
       localR2 <
         (InsideRadius + kCarTolerance) * (InsideRadius + kCarTolerance))
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
    if(std::fabs(localpos2.z()) > tubsSolid->GetZHalfLength())
      return -1;
    G4double localR2 =
      localpos2.x() * localpos2.x() + localpos2.y() * localpos2.y();
    G4double InsideRadius = tubsSolid->GetInnerRadius();
    if(localR2 >
         (InsideRadius - kCarTolerance) * (InsideRadius - kCarTolerance) &&
       localR2 <
         (InsideRadius + kCarTolerance) * (InsideRadius + kCarTolerance))
    {
      return fCurrent_Out;
    }
  }

  return -1;
}

void G4PSCylinderSurfaceCurrent::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if(HCID < 0)
    HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*) EvtMap);
}

void G4PSCylinderSurfaceCurrent::clear() { EvtMap->clear(); }

void G4PSCylinderSurfaceCurrent::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  for(const auto& [copy, current] : *(EvtMap->GetMap()))
  {
    G4cout << "  copy no.: " << copy << "  current  : ";
    if(divideByArea)
    {
      G4cout << *(current) / GetUnitValue() << " [" << GetUnit() << "]";
    }
    else
    {
      G4cout << *(current) << " [tracks]";
    }
    G4cout << G4endl;
  }
}

void G4PSCylinderSurfaceCurrent::SetUnit(const G4String& unit)
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
      G4Exception("G4PSCylinderSurfaceCurrent::SetUnit", "DetPS0002",
                  JustWarning, msg);
    }
  }
}

void G4PSCylinderSurfaceCurrent::DefineUnitAndCategory()
{
  // Per Unit Surface
  new G4UnitDefinition("percentimeter2", "percm2", "Per Unit Surface",
                       (1. / cm2));
  new G4UnitDefinition("permillimeter2", "permm2", "Per Unit Surface",
                       (1. / mm2));
  new G4UnitDefinition("permeter2", "perm2", "Per Unit Surface", (1. / m2));
}
