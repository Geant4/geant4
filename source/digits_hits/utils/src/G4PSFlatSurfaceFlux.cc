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
// $Id: G4PSFlatSurfaceFlux.cc,v 1.1 2005-11-16 23:12:42 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4PSFlatSurfaceFlux
#include "G4PSFlatSurfaceFlux.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4MultiFunctionalDetector.hh"

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

G4PSFlatSurfaceFlux::G4PSFlatSurfaceFlux(G4String name, 
					 G4int direction, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction)
{;}

G4PSFlatSurfaceFlux::~G4PSFlatSurfaceFlux()
{;}

G4bool G4PSFlatSurfaceFlux::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4VSolid * solid = 
    preStep->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
  if( solid->GetEntityType() != "G4Box" ){
    G4Exception("G4PSFlatSurfaceFluxScorer. - Solid type is not supported.");
    return FALSE;
  }
  G4Box* boxSolid = (G4Box*)(solid);

  G4int dirFlag =IsSelectedSurface(aStep,boxSolid);
  if ( dirFlag > 0 ) {
    G4int index = GetIndex(aStep);
    G4double square = 4.*boxSolid->GetXHalfLength()*boxSolid->GetXHalfLength();
    G4TouchableHandle theTouchable = preStep->GetTouchableHandle();
    G4ThreeVector pdirection = preStep->GetMomentumDirection();
    G4ThreeVector localdir  = 
      theTouchable->GetHistory()->GetTopTransform().TransformAxis(pdirection);
    G4double angleFactor = localdir.z();
    G4double flux = preStep->GetWeight(); // Current (Particle Weight)
    flux = flux/angleFactor/square;  // Flux with angle.

    if ( fDirection == fFlux_InOut || fDirection == dirFlag ){
      EvtMap->add(index,flux);
    }

#ifdef debug
    G4cout << " PASSED vol " 
	   << index << " trk "<<trkid<<" len " << fFlatSurfaceFlux<<G4endl;
#endif
  }

  return TRUE;
}

G4int G4PSFlatSurfaceFlux::IsSelectedSurface(G4Step* aStep, G4Box* boxSolid){

  G4TouchableHandle theTouchable = 
    aStep->GetPreStepPoint()->GetTouchableHandle();

  if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Entering Geometry
    G4ThreeVector stppos1= aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector localpos1 = 
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos1);
    if(fabs( localpos1.z() + boxSolid->GetZHalfLength())<kCarTolerance ){
      return fFlux_In;
    }
  }

  if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Exiting Geometry
    G4ThreeVector stppos2= aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector localpos2 = 
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos2);
    if(fabs( localpos2.z() + boxSolid->GetZHalfLength())<kCarTolerance ){
      return fFlux_Out;
    }
  }

  return -1;
}

void G4PSFlatSurfaceFlux::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
				    GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSFlatSurfaceFlux::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSFlatSurfaceFlux::clear(){
  EvtMap->clear();
}

void G4PSFlatSurfaceFlux::DrawAll()
{;}

void G4PSFlatSurfaceFlux::PrintAll()
{
  G4cout << " PrimitiveSenstivity " << GetName() <<G4endl; 
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  flux  : " << *(itr->second) /mm2
	   << G4endl;
  }
}

