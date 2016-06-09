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
// $Id: G4PSSphereSurfaceCurrent.cc,v 1.2 2007/08/14 21:23:52 taso Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// G4PSSphereSurfaceCurrent
#include "G4PSSphereSurfaceCurrent.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4GeometryTolerance.hh"
////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring only Surface Current.
//  Current version assumes only for G4Sphere shape. 
//
// Surface is defined  at the inside of sphere.
// Direction                  -Rmin   +Rmax
//   0  IN || OUT            ->|<-     |
//   1  IN                   ->|       |
//   2  OUT                    |<-     |
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 
///////////////////////////////////////////////////////////////////////////////

G4PSSphereSurfaceCurrent::G4PSSphereSurfaceCurrent(G4String name, 
					 G4int direction, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction),
     weighted(true),divideByArea(true)
{;}

G4PSSphereSurfaceCurrent::~G4PSSphereSurfaceCurrent()
{;}

G4bool G4PSSphereSurfaceCurrent::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4VSolid * solid = 
    preStep->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
  if( solid->GetEntityType() != "G4Sphere" ){
    G4Exception("G4PSSphereSurfaceCurrentScorer. - Solid type is not supported.");
    return FALSE;
  }
  G4Sphere* sphereSolid = (G4Sphere*)(solid);

  G4int dirFlag =IsSelectedSurface(aStep,sphereSolid);
  if ( dirFlag > 0 ) {
    if ( fDirection == fCurrent_InOut || fDirection == dirFlag ){
	G4double radi   = sphereSolid->GetInsideRadius();
	G4double dph    = sphereSolid->GetDeltaPhiAngle()/radian;
	G4double stth   = sphereSolid->GetStartThetaAngle()/radian;
	G4double enth   = stth+sphereSolid->GetDeltaThetaAngle()/radian;
	G4double current = 1.0;
	if ( weighted) current = preStep->GetWeight(); // Current (Particle Weight)
	if ( divideByArea ){
	    G4double square = radi*radi*dph*( -std::cos(enth) + std::cos(stth) );
	    current = current/square;  // Current with angle.
	}

	G4int index = GetIndex(aStep);
	EvtMap->add(index,current);
    }
  }

  return TRUE;
}

G4int G4PSSphereSurfaceCurrent::IsSelectedSurface(G4Step* aStep, G4Sphere* sphereSolid){

  G4TouchableHandle theTouchable = 
    aStep->GetPreStepPoint()->GetTouchableHandle();
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  
  if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Entering Geometry
    G4ThreeVector stppos1= aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector localpos1 = 
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos1);
    G4double localR2 = localpos1.x()*localpos1.x()
                      +localpos1.y()*localpos1.y()
                      +localpos1.z()*localpos1.z();
    G4double InsideRadius2 = 
      sphereSolid->GetInsideRadius()*sphereSolid->GetInsideRadius();
    if(std::fabs( localR2 - InsideRadius2 ) < kCarTolerance ){
      return fCurrent_In;
    }
  }

  if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Exiting Geometry
    G4ThreeVector stppos2= aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector localpos2 = 
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos2);
    G4double localR2 = localpos2.x()*localpos2.x()
                      +localpos2.y()*localpos2.y()
                      +localpos2.z()*localpos2.z();
    G4double InsideRadius2 = 
      sphereSolid->GetInsideRadius()*sphereSolid->GetInsideRadius();
    if(std::fabs( localR2 - InsideRadius2 ) < kCarTolerance ){
      return fCurrent_Out;
    }
  }

  return -1;
}

void G4PSSphereSurfaceCurrent::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSSphereSurfaceCurrent::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSSphereSurfaceCurrent::clear(){
  EvtMap->clear();
}

void G4PSSphereSurfaceCurrent::DrawAll()
{;}

void G4PSSphereSurfaceCurrent::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() <<G4endl; 
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first  << "  current  : " ;
    if ( divideByArea ) G4cout << *(itr->second)*cm*cm << " [/cm2]" ;
    else G4cout << *(itr->second)*cm*cm << " [track]" ;
    G4cout  << G4endl;
  }
}

