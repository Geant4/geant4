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

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SemiInfiniteTarget.hh"
#include "AddOnTargetLayer.hh"
#include "RegionInformation.hh"
#include "TargetGeometryManager.hh"
#include "AnalysisManager.hh"
#include "Materials.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RegionStore.hh"
#include "G4VisAttributes.hh"


DetectorConstruction::DetectorConstruction() : 
    worldVolPhys(0) {

  updateManager = new TargetGeometryManager();
  messenger = new DetectorMessenger(this);

  targetRegionName = "Target";
  calThickness = 1.0 * mm;
}


DetectorConstruction::~DetectorConstruction() {
 
  delete messenger;
  delete updateManager;
}


G4VPhysicalVolume* DetectorConstruction::Construct() {

  G4Material* worldVolMat = 
      Materials::Instance() -> GetMaterial("Galactic");   

  G4Box* worldVolSolid = new G4Box("World",50 * cm, 50 * cm, 50 * cm);
  
  G4LogicalVolume* worldVolLogic = new G4LogicalVolume(worldVolSolid, 
                                                       worldVolMat, "World");
  worldVolLogic -> SetVisAttributes (G4VisAttributes::Invisible);

  worldVolPhys = new G4PVPlacement(0, G4ThreeVector(), "World", 
                                   worldVolLogic, 0, false, 0);

  RegionInformation* regionInfoWorld = new RegionInformation();
  regionInfoWorld -> FlagRegionAsWorld();

  G4Region* worldRegion = (*(G4RegionStore::GetInstance()))[0];
  worldRegion -> SetUserInformation(regionInfoWorld); 


  SemiInfiniteTarget* infiniteSlab = 
      new SemiInfiniteTarget("infiniteTargetLayer", updateManager, 
                             worldVolPhys);

  RegionInformation* regionInfoTarget = new RegionInformation();
  regionInfoTarget -> FlagRegionAsTarget();

  G4Region* targetRegion = new G4Region(targetRegionName);
  targetRegion -> SetUserInformation(regionInfoTarget);

  targetRegion -> AddRootLogicalVolume(infiniteSlab -> GetVolLog());

  return worldVolPhys;
}


void DetectorConstruction::CreateFrontLayer(G4String layerName) {

  TargetComponent* currentFrontLayer = updateManager -> GetFrontLayer();

  AddOnTargetLayer* frontLayer = 
      new AddOnTargetLayer(currentFrontLayer, updateManager, layerName);

  G4Region* targetRegion = 
      G4RegionStore::GetInstance() -> GetRegion(targetRegionName);
  targetRegion -> AddRootLogicalVolume(frontLayer -> GetVolLog());  

  G4RunManager::GetRunManager() -> DefineWorldVolume(worldVolPhys);
}


void DetectorConstruction::SetLayerRadius(G4double rad) {

  updateManager -> SetRadius(rad);
  G4RunManager::GetRunManager() -> DefineWorldVolume(worldVolPhys);
}


void DetectorConstruction::SetLayerThickness(G4double thickn) {

  updateManager -> SetThickness(thickn);
  G4RunManager::GetRunManager() -> DefineWorldVolume(worldVolPhys);
}


void DetectorConstruction::SetLayerMaterial(G4String mat) {

  updateManager -> SetMaterial(mat);
  G4RunManager::GetRunManager() -> DefineWorldVolume(worldVolPhys);
}


void DetectorConstruction::SetLayerMaxStepSize(G4double max) {

  updateManager -> SetMaxStepSize(max);
  G4RunManager::GetRunManager() -> DefineWorldVolume(worldVolPhys);
}


void DetectorConstruction::CreateCalorimeter(G4double zPosition) {

  AnalysisManager::Instance() -> 
      CreateCalorimeter(zPosition, calThickness,updateManager -> GetRadius());
}


void DetectorConstruction::SetCalorimeterThickness(G4double thickn) {

  if(thickn > 0.0 * mm) {
     calThickness = thickn;
  }
}
