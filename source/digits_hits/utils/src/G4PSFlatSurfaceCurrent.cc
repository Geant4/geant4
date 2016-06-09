//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PSFlatSurfaceCurrent.cc,v 1.3 2005/11/19 03:16:07 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// G4PSFlatSurfaceCurrent
#include "G4PSFlatSurfaceCurrent.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
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
// 
///////////////////////////////////////////////////////////////////////////////


G4PSFlatSurfaceCurrent::G4PSFlatSurfaceCurrent(G4String name, 
					 G4int direction, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction)
{;}

G4PSFlatSurfaceCurrent::~G4PSFlatSurfaceCurrent()
{;}

G4bool G4PSFlatSurfaceCurrent::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4VSolid * solid = 
    preStep->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
  if( solid->GetEntityType() != "G4Box" ){
    G4Exception("G4PSFlatSurfaceCurrentScorer. - Solid type is not supported.");
    return FALSE;
  }
  G4Box* boxSolid = (G4Box*)(solid);

  G4int dirFlag =IsSelectedSurface(aStep,boxSolid);
  if ( dirFlag > 0 ) {
    G4int index = GetIndex(aStep);
    G4double square = 4.*boxSolid->GetXHalfLength()*boxSolid->GetYHalfLength();
    G4TouchableHandle theTouchable = preStep->GetTouchableHandle();
    G4ThreeVector pdirection = preStep->GetMomentumDirection();
    G4ThreeVector localdir  = 
      theTouchable->GetHistory()->GetTopTransform().TransformAxis(pdirection);
    G4double current = preStep->GetWeight(); // Current (Particle Weight)
    current = current/square;  // Current with angle.

    if ( fDirection == fCurrent_InOut || fDirection == dirFlag ){
      EvtMap->add(index,current);
    }

  }

  return TRUE;
}

G4int G4PSFlatSurfaceCurrent::IsSelectedSurface(G4Step* aStep, G4Box* boxSolid){

  G4TouchableHandle theTouchable = 
    aStep->GetPreStepPoint()->GetTouchableHandle();
  

  if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Entering Geometry
    G4ThreeVector stppos1= aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector localpos1 = 
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos1);
    if(fabs( localpos1.z() + boxSolid->GetZHalfLength())<kCarTolerance ){
      return fCurrent_In;
    }
  }

  if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Exiting Geometry
    G4ThreeVector stppos2= aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector localpos2 = 
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos2);
    if(fabs( localpos2.z() + boxSolid->GetZHalfLength())<kCarTolerance ){
      return fCurrent_Out;
    }
  }

  return -1;
}

void G4PSFlatSurfaceCurrent::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSFlatSurfaceCurrent::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSFlatSurfaceCurrent::clear(){
  EvtMap->clear();
}

void G4PSFlatSurfaceCurrent::DrawAll()
{;}

void G4PSFlatSurfaceCurrent::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() <<G4endl; 
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  current  : " << *(itr->second)*cm*cm << " [cm^-2]"
	   << G4endl;
  }
}

