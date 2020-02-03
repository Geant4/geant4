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
// // G4PSCylinderSurfaceFlux
#include "G4PSCylinderSurfaceFlux.hh"

#include "G4SystemOfUnits.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
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
// 2010-07-22   Introduce Unit specification.
// 2010-07-22   Add weighted and divideByArea options
// 2011-02-21   Get correct momentum direction in Flux_Out.
///////////////////////////////////////////////////////////////////////////////

G4PSCylinderSurfaceFlux::G4PSCylinderSurfaceFlux(G4String name, 
						 G4int direction, G4int depth)
    : G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction),EvtMap(0),
      weighted(true),divideByArea(true)
{
    DefineUnitAndCategory();
    SetUnit("percm2");
}

G4PSCylinderSurfaceFlux::G4PSCylinderSurfaceFlux(G4String name, 
						 G4int direction, 
						 const G4String& unit, 
						 G4int depth)
    : G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction),EvtMap(0),
      weighted(true),divideByArea(true)
{
    DefineUnitAndCategory();
    SetUnit(unit);
}

G4PSCylinderSurfaceFlux::~G4PSCylinderSurfaceFlux()
{;}

G4bool G4PSCylinderSurfaceFlux::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4StepPoint* preStep = aStep->GetPreStepPoint();

  G4VPhysicalVolume* physVol = preStep->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid * solid = 0;
  if(physParam)
  { // for parameterized volume
    G4int idx = ((G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable()))
                ->GetReplicaNumber(indexDepth);
    solid = physParam->ComputeSolid(idx, physVol);
    solid->ComputeDimensions(physParam,idx,physVol);
  }
  else
  { // for ordinary volume
    solid = physVol->GetLogicalVolume()->GetSolid();
  }

  G4Tubs* tubsSolid = (G4Tubs*)(solid);
  
  G4int dirFlag =IsSelectedSurface(aStep,tubsSolid);
  
  if ( dirFlag > 0 ){
    if (fDirection == fFlux_InOut || dirFlag == fDirection ){

      G4StepPoint* thisStep=0;
      if ( dirFlag == fFlux_In ){
	thisStep = preStep;
      }else if ( dirFlag == fFlux_Out ){
	thisStep = aStep->GetPostStepPoint();
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
    
      G4double flux = 1.0;
      if ( weighted ) flux *=preStep->GetWeight();  
      // Current (Particle Weight)

      flux = flux/angleFactor;   
      if ( divideByArea ) flux /= square;
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
    //if(std::fabs( localpos1.x()*localpos1.x()+localpos1.y()*localpos1.y()  
    //   - (tubsSolid->GetInnerRadius()*tubsSolid->GetInnerRadius()))
    //   <kCarTolerance ){
    G4double localR2 = localpos1.x()*localpos1.x()+localpos1.y()*localpos1.y();
    G4double InsideRadius = tubsSolid->GetInnerRadius();
    if (localR2 > (InsideRadius-kCarTolerance)*(InsideRadius-kCarTolerance)
	&&localR2 < (InsideRadius+kCarTolerance)*(InsideRadius+kCarTolerance)){
      return fFlux_In;
    }
  }

  if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Exiting Geometry
    G4ThreeVector stppos2= aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector localpos2 = 
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos2);
    if ( std::fabs(localpos2.z()) > tubsSolid->GetZHalfLength() ) return -1;
    //if(std::fabs( localpos2.x()*localpos2.x()+localpos2.y()*localpos2.y()  
    //   - (tubsSolid->GetInnerRadius()*tubsSolid->GetInnerRadius())) 
    // <kCarTolerance ){
    G4double localR2 = localpos2.x()*localpos2.x()+localpos2.y()*localpos2.y();
    G4double InsideRadius = tubsSolid->GetInnerRadius();
    if (localR2 > (InsideRadius-kCarTolerance)*(InsideRadius-kCarTolerance)
	&&localR2 < (InsideRadius+kCarTolerance)*(InsideRadius+kCarTolerance)){
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
	   << "  flux  : " << *(itr->second)/GetUnitValue()
	   << " ["<<GetUnit()<<"]"
	   << G4endl;
  }
}

void G4PSCylinderSurfaceFlux::SetUnit(const G4String& unit)
{
    if ( divideByArea ) {
	CheckAndSetUnit(unit,"Per Unit Surface");
    } else {
	if (unit == "" ){
	    unitName = unit;
	    unitValue = 1.0;
	}else{
	    G4String msg = "Invalid unit ["+unit+"] (Current  unit is [" +GetUnit()+"] ) for " + GetName();
	    G4Exception("G4PSCylinderSurfaceFlux::SetUnit","DetPS0003",JustWarning,msg);
	}
    }
}

void G4PSCylinderSurfaceFlux::DefineUnitAndCategory(){
   // Per Unit Surface
   new G4UnitDefinition("percentimeter2","percm2","Per Unit Surface",(1./cm2));
   new G4UnitDefinition("permillimeter2","permm2","Per Unit Surface",(1./mm2));
   new G4UnitDefinition("permeter2","perm2","Per Unit Surface",(1./m2));
}


