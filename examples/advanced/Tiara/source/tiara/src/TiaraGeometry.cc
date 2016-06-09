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
// $Id: TiaraGeometry.cc,v 1.3 2003/06/25 09:13:03 gunter Exp $
// GEANT4 tag $Name: geant4-06-00 $
//

#include "TiaraGeometry.hh"
#include "TiaraDimensions.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "TiaraMaterials.hh"
#include "G4PVPlacement.hh"
#include "TiaraConcreteShieldA.hh"
#include "TiaraIronShieldA.hh"
#include "TiaraConcreteShieldB.hh"
#include "TiaraIronShieldB.hh"
#include "TiaraPart.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

TiaraGeometry::TiaraGeometry(TiaraMaterials &matfac) :
  fLogicalWorld(0),
  fWorldVolume(0),
  fMaterials(matfac),
  fShieldWidth(0),
  fColimatorWidth(0),
  fColliPipeLogVol(0),
  fTiaraDimensions(TiaraDimensions())
{}

TiaraGeometry::~TiaraGeometry()
{}
  
G4VPhysicalVolume* TiaraGeometry::Construct(){
  PlaceComponents();
  return fWorldVolume;
}

G4Material *TiaraGeometry::GetWorldMaterial() {
  return fMaterials.GetMaterial("vacuum");
}


void TiaraGeometry::ConstHall(){
  G4VSolid *worldSolid = new G4Box("worldBox", 
				   fTiaraDimensions.worldHalfWidth,
				   fTiaraDimensions.worldHalfWidth,
				   fTiaraDimensions.worldHalfLength);
  fLogicalWorld = 
    new G4LogicalVolume(worldSolid, GetWorldMaterial(), "LogicalWorld");
  G4String name("Hall");
  fWorldVolume = new 
    G4PVPlacement(0, G4ThreeVector(0,0,0), fLogicalWorld,
		  name, 0, false, 0);
}

void TiaraGeometry::CreateComponents(){

  fTiaraComponents.
    insert(new TiaraConcreteShieldA(fMaterials, fTiaraDimensions));

  fTiaraComponents.
    insert(new TiaraIronShieldA(fMaterials, fTiaraDimensions));

  fTiaraComponents.
    insert(new TiaraConcreteShieldB(fMaterials, fTiaraDimensions));
  
  fTiaraComponents.
    insert(new TiaraIronShieldB(fMaterials, fTiaraDimensions));

}

G4VPhysicalVolume *TiaraGeometry::
AddPhysicalDetector(G4double xDist,
		    const G4String &physName) {
  G4Tubs *tub = new G4Tubs("detectorTub",
			   0,
			   fTiaraDimensions.detectorRadius,
			   fTiaraDimensions.detectorHalfHight,
			   0,
			   360 * deg);
  
  G4ThreeVector pos(xDist, 0,  fTiaraDimensions.targetPosZ + 
		    fTiaraDimensions.distTargetExperiment + 
		    fShieldWidth + fColimatorWidth + 
		    fTiaraDimensions.detectorHalfHight);
  
  G4Material *mat = fMaterials.GetMaterial("vacuum");
  G4LogicalVolume *logVol = 
    new G4LogicalVolume(tub, mat, "logDetectorTub");
  G4VisAttributes* pVis = new 
    G4VisAttributes(G4Colour(1, 0, 0));
  logVol->SetVisAttributes(pVis);
  

  return PlaceExpComponent(pos, logVol, physName);

}

G4VPhysicalVolume *TiaraGeometry::
AddPhysicalRingDetector(G4double xDist,
			const G4String &physName){
  if (xDist < fTiaraDimensions.detectorRadius) {
    G4cout << "TiaraGeometry::AddPhysicalRingDetector: (xDist < fTiaraDimensions.detectorRadius" << G4endl;
    return 0;
  }
  G4Tubs *tub = new G4Tubs("detectorTub",
			   xDist - fTiaraDimensions.detectorRadius,
			   xDist + fTiaraDimensions.detectorRadius,
			   fTiaraDimensions.detectorHalfHight,
			   0,
			   360 * deg);
  G4ThreeVector pos(0, 0,  fTiaraDimensions.targetPosZ + 
		    fTiaraDimensions.distTargetExperiment + 
		    fShieldWidth + fColimatorWidth + 
		    fTiaraDimensions.detectorHalfHight);

  G4Material *mat = fMaterials.GetMaterial("vacuum");
  G4LogicalVolume *logVol = 
    new G4LogicalVolume(tub, mat, "logRingDetectorTub");
  G4VisAttributes* pVis = new 
    G4VisAttributes(G4Colour(1, 0, 0));
  logVol->SetVisAttributes(pVis);

  return PlaceExpComponent(pos, logVol, physName);
  
}

G4VPhysicalVolume *TiaraGeometry::
AddDetectorSlab( const G4String &physName) {
  G4Tubs *tub = new G4Tubs("detectorTub",
			   0,
			   0.5*fTiaraDimensions.widthExperiment,
			   fTiaraDimensions.detectorHalfHight,
			   0,
			   360 * deg);
  G4ThreeVector pos(0, 0,  fTiaraDimensions.targetPosZ + 
		    fTiaraDimensions.distTargetExperiment + 
		    fShieldWidth + fColimatorWidth + 
		    fTiaraDimensions.detectorHalfHight);

  G4Material *mat = fMaterials.GetMaterial("vacuum");
  G4LogicalVolume *logVol = 
    new G4LogicalVolume(tub, mat, "logSlabDetectorTub");
  G4VisAttributes* pVis = new 
    G4VisAttributes(G4Colour(1, 0, 0));
  logVol->SetVisAttributes(pVis);

  return PlaceExpComponent(pos, logVol, physName);
  
}

G4VPhysicalVolume *TiaraGeometry::
AddSourceDetector() {
  G4Tubs *tub = new G4Tubs("srcTub",
			   0,
			   fTiaraDimensions.pipeRadius,
			   fTiaraDimensions.srcDetectorWidth / 2,
			   0,
			   360 * deg);
  G4Material *mat = fMaterials.GetMaterial("vacuum");
  G4LogicalVolume *logVol = 
    new G4LogicalVolume(tub, mat, "logSrcDetector");

  G4VisAttributes* pVis = new 
    G4VisAttributes(G4Colour(1, 0, 0));
  logVol->SetVisAttributes(pVis);
  
  G4VPhysicalVolume *physVol = 0;
  if (fColliPipeLogVol) {
    G4ThreeVector pos(0, 0,  fColimatorWidth/2 - 
		      fTiaraDimensions.srcDetectorWidth / 2);
    physVol = new G4PVPlacement(0,  
				pos, 
				logVol,
				"sourceDetector", 
				fColliPipeLogVol, 
				false, 
				0);
  }
  else {
    G4ThreeVector pos(0, 0,  fTiaraDimensions.targetPosZ + 
		      fTiaraDimensions.distTargetExperiment + 
		      fColimatorWidth - 
		      fTiaraDimensions.srcDetectorWidth / 2);
    physVol = PlaceExpComponent(pos, logVol, "sourceDetector");
  }

  return physVol;
}


void TiaraGeometry::PlaceComponents(){
  for (TiaraComponents::iterator it = fTiaraComponents.begin();
       it != fTiaraComponents.end(); ++it) {
    TiaraParts tiaraParts = (*it)->GetParts();
    for (TiaraParts::iterator ip = tiaraParts.begin();
	 ip != tiaraParts.end(); ++ip){
      TiaraPart tiaraPart = *ip;
      G4LogicalVolume *logVol = ip->logVol;
      G4cout << logVol->GetName() 
	     << ", pos: " << ip->pos << G4endl;
      G4String name(logVol->GetName());
      name += "_placed";
      new G4PVPlacement(ip->rot,  
			ip->pos, 
			logVol,
			name, 
			fLogicalWorld, 
			false, 
			0);
    }
  }
}


void TiaraGeometry::BuildGeometry(const TiaraDimensions &td){
  fTiaraDimensions = td;
  ConstHall();
}

G4LogicalVolume *TiaraGeometry::BuildShield(G4double width,
					    const G4String &matName){

  fShieldWidth = width;

  G4LogicalVolume *logVol = 0;
  if (! (width > 0)) {
    G4cout << "TiaraGeometry::BuildShiel witdth must be > 0"
	   << G4endl;
  }
  else {
    G4Material *mat = 0;
    mat = fMaterials.GetMaterial(matName);
    if (!mat) {
      G4cout << "TiaraGeometry::BuildShiel no valid material"
	     << G4endl;
    }
    else {
      G4Box  *box = new G4Box("box",
			      fTiaraDimensions.widthExperiment/2,
			      fTiaraDimensions.widthExperiment/2,
			      width/2);
      logVol = new G4LogicalVolume(box, mat, "logShieldBox");
      G4VisAttributes* pVis = new 
	G4VisAttributes(G4Colour(1, 1, 0));
      logVol->SetVisAttributes(pVis);
    }
  }
  return logVol;
}

G4LogicalVolume *TiaraGeometry::
BuildCollimator(G4double width,
		const G4String &outerMatName,
		const G4String &innerMatName) {
  fColimatorWidth = width;
    
  G4LogicalVolume *logVol = 0;
  if (! (width > 0)) {
    G4cout << "TiaraGeometry::BuildCollimator witdth must be > 0"
	   << G4endl;
  }
  else {
    G4Material *matOuter = 0;
    matOuter = fMaterials.GetMaterial(outerMatName);
    G4Material *matInner = 0;
    matInner = fMaterials.GetMaterial(innerMatName);
    if (!matOuter || !matInner) {
      G4cout << "TiaraGeometry::BuildCollimator no valid material"
	     << G4endl;
    }
    else {
      G4Box  *box = new G4Box("box",
			      fTiaraDimensions.widthExperiment/2,
			      fTiaraDimensions.widthExperiment/2,
			      width/2);
      logVol = new G4LogicalVolume(box, matOuter, "logColimatorBox");
      G4VisAttributes* pVis = new 
	G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
      logVol->SetVisAttributes(pVis);
      
      G4Tubs *tub = new G4Tubs("tub",
			       0,
			       fTiaraDimensions.pipeRadius,
			       width/2,
			       0,
			       360 * deg);
      
      fColliPipeLogVol = new G4LogicalVolume(tub, matInner,
					    "logColimatorTub");
      
      new G4PVPlacement(0,  
			G4ThreeVector(0, 0, 0),
			fColliPipeLogVol,
			"pipeInCollimator", 
			logVol, 
			false, 
			0);
    }
  }
  return logVol;
}

G4VPhysicalVolume *TiaraGeometry::
PlaceExpComponent(const G4ThreeVector &pos, 
		  G4LogicalVolume *logVol,
		  const G4String &physName){
  
  return new G4PVPlacement(0,  
			   pos, 
			   logVol,
			   physName, 
			   fLogicalWorld, 
			   false, 
			   0);
  
}

