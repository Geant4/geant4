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
// $Id: G4PSCylinderSurfaceFlux.cc,v 1.2 2008/12/18 12:57:18 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// // G4PSCylinderSurfaceFlux
#include "G4PSCylinderSurfaceFlux.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4GeometryTolerance.hh"
// ////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring Surface Flux.
//  Current version assumes only for G4Tubs shape, and the surface
//  is fixed on inner plane of the tube.
//
// Surface is defined at the innner surface of the tube.
// Direction                   R    R+dR
//   0  IN || OUT            ->|<-  |
//   1  IN                   ->|    |
//   2  OUT                    |<-  |
//
// Created: 2007-03-29  Tsukasa ASO
///////////////////////////////////////////////////////////////////////////////

G4PSCylinderSurfaceFlux::G4PSCylinderSurfaceFlux(G4String name, 
						 G4int direction, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction)
{;}

G4PSCylinderSurfaceFlux::~G4PSCylinderSurfaceFlux()
{;}

G4bool G4PSCylinderSurfaceFlux::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4StepPoint* postStep = aStep->GetPreStepPoint();
  G4VSolid * solid = 
    preStep->GetPhysicalVolume()->GetLogicalVolume()->GetSolid();
  if( solid->GetEntityType() != "G4Tubs" ){
    G4Exception("G4PSCylinderSurfaceFluxScorer. - Solid type is not supported.");
    return FALSE;
  }
  G4Tubs* tubsSolid = (G4Tubs*)(solid);
  
  G4int dirFlag =IsSelectedSurface(aStep,tubsSolid);
  
  if ( dirFlag > 0 ){
    if (fDirection == fFlux_InOut || dirFlag == fDirection ){

      G4StepPoint* thisStep=0;
      if ( dirFlag == fFlux_In ){
	thisStep = preStep;
      }else if ( dirFlag == fFlux_Out ){
	thisStep = postStep;
      }else{
	return FALSE;
      }
  
      G4TouchableHandle theTouchable = thisStep->GetTouchableHandle();
      G4ThreeVector pdirection = thisStep->GetMomentumDirection();
      G4ThreeVector localdir  = 
	theTouchable->GetHistory()->GetTopTransform().TransformAxis(pdirection);
      G4ThreeVector position = thisStep->GetPosition();
      G4ThreeVector localpos  =
	theTouchable->GetHistory()->GetTopTransform().TransformAxis(position);
      G4double angleFactor = (localdir.x()*localpos.x()+localdir.y()*localpos.y())
	/std::sqrt(localdir.x()*localdir.x()
		   +localdir.y()*localdir.y()+localdir.z()*localdir.z())
	/std::sqrt(localpos.x()*localpos.x()+localpos.y()*localpos.y());
    
      if ( angleFactor < 0 ) angleFactor *= -1.;
      G4double square = 2.*tubsSolid->GetZHalfLength()
	*tubsSolid->GetInnerRadius()* tubsSolid->GetDeltaPhiAngle()/radian;
    
      G4double flux = preStep->GetWeight();  
      // Current (Particle Weight)

      flux = flux/angleFactor/square;   
      //Flux with angle.
      G4int index = GetIndex(aStep);
      EvtMap->add(index,flux);
      return TRUE;
    }else{
      return FALSE;
    }
  }else{
      return FALSE;
  }
}

G4int G4PSCylinderSurfaceFlux::IsSelectedSurface(G4Step* aStep, G4Tubs* tubsSolid){

  G4TouchableHandle theTouchable = 
    aStep->GetPreStepPoint()->GetTouchableHandle();
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();  

  if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Entering Geometry
    G4ThreeVector stppos1= aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector localpos1 = 
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos1);
    if ( std::fabs(localpos1.z()) > tubsSolid->GetZHalfLength() ) return -1;
    if(std::fabs( localpos1.x()*localpos1.x()+localpos1.y()*localpos1.y() ) 
       - (tubsSolid->GetInnerRadius()*tubsSolid->GetInnerRadius())<kCarTolerance ){
      return fFlux_In;
    }
  }

  if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Exiting Geometry
    G4ThreeVector stppos2= aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector localpos2 = 
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos2);
    if ( std::fabs(localpos2.z()) > tubsSolid->GetZHalfLength() ) return -1;
    if(std::fabs( localpos2.x()*localpos2.x()+localpos2.y()*localpos2.y() ) 
       - (tubsSolid->GetInnerRadius()*tubsSolid->GetInnerRadius())<kCarTolerance ){
      return fFlux_Out;
    }
  }

  return -1;
}

void G4PSCylinderSurfaceFlux::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
				    GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSCylinderSurfaceFlux::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSCylinderSurfaceFlux::clear(){
  EvtMap->clear();
}

void G4PSCylinderSurfaceFlux::DrawAll()
{;}

void G4PSCylinderSurfaceFlux::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer" << GetName() <<G4endl; 
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  flux  : " << *(itr->second)*cm*cm << " [cm^-2]"
	   << G4endl;
  }
}

