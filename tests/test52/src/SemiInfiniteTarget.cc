
#include "SemiInfiniteTarget.hh"
#include "TargetGeometryManager.hh"
#include "TargetLayerSD.hh"
#include "Materials.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"


SemiInfiniteTarget::SemiInfiniteTarget(G4String layerName,
                                       TargetGeometryManager* geomManager,
                                       G4VPhysicalVolume* world,
                                       G4String material,
                                       G4double thickn,
                                       G4double rad,
                                       G4double maxStep) :
    TargetComponent(geomManager, layerName, thickn, material),
    worldVolPhys(world), 
    radius(rad),
    maxStepSize(maxStep),
    semiInfTargetSolid(0),
    semiInfTargetLogic(0),
    semiInfTargetPhys(0) {

  updateManager -> Attach(this);  

  semiInfTargetSolid = new G4Tubs(Name(), 0.0 * mm, radius, Thickness() * 0.5,
				  0.0 * deg, 360.0 * deg); 
      
  semiInfTargetLogic = new G4LogicalVolume(semiInfTargetSolid,
                                           Material(),
                                           Name());

  semiInfTargetPhys = new G4PVPlacement(0, 
                         G4ThreeVector(0.0 * cm, 0.0 * cm, Thickness() * 0.5),
                         Name(), semiInfTargetLogic, worldVolPhys, false, 0);

  TargetLayerSD* semiInfTargetSD = new TargetLayerSD(Name()); 

  semiInfTargetLogic -> SetSensitiveDetector(semiInfTargetSD);

  G4SDManager* SDManager = G4SDManager::GetSDMpointer();
  SDManager -> AddNewDetector(semiInfTargetSD);

  semiInfTargetUserLimits = new G4UserLimits(maxStepSize);
  semiInfTargetLogic -> SetUserLimits(semiInfTargetUserLimits);

  G4VisAttributes* semiInfTargetVisAtt = 
                      new G4VisAttributes(G4Colour::Yellow());
  semiInfTargetVisAtt -> SetVisibility(true);
  semiInfTargetLogic -> SetVisAttributes(semiInfTargetVisAtt);

  GeometryUpdate(updateManager);
}


SemiInfiniteTarget::~SemiInfiniteTarget() {

}


void SemiInfiniteTarget::GeometryUpdate(TargetGeometryManager* geomManager) {

  if(geomManager == updateManager) {

     G4double newRadius = updateManager -> GetRadius();
     if(newRadius > 0.0 * cm) {
	radius = newRadius;
        semiInfTargetSolid -> SetOuterRadius(radius);
     } 

     G4double newMaxStepSize = updateManager -> GetMaxStepSize();
     if(newMaxStepSize > 0.0 * mm) {
        maxStepSize = newMaxStepSize;
        semiInfTargetUserLimits -> SetMaxAllowedStep(maxStepSize);
     }

     Thickness(updateManager -> GetThickness(this));
     semiInfTargetSolid -> SetZHalfLength(Thickness() * 0.5);

     Material(updateManager -> GetMaterial(this));
     semiInfTargetLogic -> SetMaterial(Material());     

     G4double position = updateManager -> GetPosition(this);
     semiInfTargetPhys -> SetTranslation(
         G4ThreeVector(0.0 * cm, 0.0 * cm, position));
  }
}


