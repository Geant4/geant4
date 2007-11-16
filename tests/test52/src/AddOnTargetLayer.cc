
#include "AddOnTargetLayer.hh"
#include "TargetGeometryManager.hh"
#include "TargetLayerSD.hh"
#include "Materials.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"


AddOnTargetLayer::AddOnTargetLayer(TargetComponent* component,
                                   TargetGeometryManager* geomManager,
                                   G4String layerName,
                                   G4String material,
                                   G4double thickn) : 
    TargetComponent(geomManager, layerName, thickn, material),
    targetComponent(component) {

  updateManager -> Attach(this);

  addOnTargetLayerSolid = new G4Tubs(Name(),
                                  0.0 * mm, targetComponent -> GetRadius(), 
				  Thickness() * 0.5, 0.0 * deg, 360.0 * deg); 
      
  addOnTargetLayerLogic = new G4LogicalVolume(
                                           addOnTargetLayerSolid,
                                           Material(),
                                           Name());

  addOnTargetLayerPhys = new G4PVPlacement(0, 
                          G4ThreeVector(0.0 * cm, 0.0 * cm, Thickness() * 0.5),
                          Name(), addOnTargetLayerLogic,
                          targetComponent -> GetWorldVolPhys(), false, 0);

  TargetLayerSD* addOnTargetLayerSD = 
                      new TargetLayerSD(Name()); 
  addOnTargetLayerLogic -> SetSensitiveDetector(addOnTargetLayerSD);

  G4SDManager* SDManager = G4SDManager::GetSDMpointer();
  SDManager -> AddNewDetector(addOnTargetLayerSD);

  addOnTargetLayerUserLimits = 
                      new G4UserLimits(targetComponent -> GetMaxStepSize());
  addOnTargetLayerLogic -> SetUserLimits(addOnTargetLayerUserLimits);

  G4VisAttributes* addOnTargetLayerVisAtt = 
                      new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
  addOnTargetLayerVisAtt -> SetVisibility(true);
  addOnTargetLayerLogic -> SetVisAttributes(addOnTargetLayerVisAtt);
}


AddOnTargetLayer::~AddOnTargetLayer() {

}


void AddOnTargetLayer::GeometryUpdate(TargetGeometryManager* geomManager) {
 
  if(geomManager == updateManager) {

     G4double newRadius = updateManager -> GetRadius();
     if(newRadius > 0.0 * cm) {
        addOnTargetLayerSolid -> SetOuterRadius(newRadius);
     } 

     G4double newMaxStepSize = updateManager -> GetMaxStepSize();
     if(newMaxStepSize > 0.0 * mm) {
        addOnTargetLayerUserLimits -> SetMaxAllowedStep(newMaxStepSize);
     }
 
     Thickness(updateManager -> GetThickness(this));
     addOnTargetLayerSolid -> SetZHalfLength(Thickness() * 0.5);

     Material(updateManager -> GetMaterial(this));
     addOnTargetLayerLogic -> SetMaterial(Material());     

     G4double position = updateManager -> GetPosition(this);
     addOnTargetLayerPhys -> SetTranslation(
        G4ThreeVector(0.0 * cm, 0.0 * cm, position));
  }
}
