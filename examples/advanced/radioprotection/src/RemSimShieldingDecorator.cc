#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "RemSimShieldingDecorator.hh"
#include "RemSimDecorator.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "RemSimSensitiveDetector.hh"
#include "RemSimROGeometry.hh"
#include "G4SDManager.hh"

RemSimShieldingDecorator::
RemSimShieldingDecorator(RemSimVGeometryComponent* comp)
  : RemSimDecorator(comp)
{
  
}
RemSimShieldingDecorator::~RemSimShieldingDecorator()
{
  ;
}
void RemSimShieldingDecorator::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
   RemSimDecorator::ConstructComponent(motherVolume);
  ConstructShielding(motherVolume);
}

void RemSimShieldingDecorator::DestroyComponent()
{
  ;
}
void RemSimShieldingDecorator::ConstructShielding(G4VPhysicalVolume* mother)
{
  // Materials for the astronaut composition
  RemSimMaterial* pMaterial = new RemSimMaterial();
  G4Material* water = pMaterial -> GetMaterial("Water");
  delete pMaterial;

  // Geometry definition
  G4double shieldingX = 5.*m;
  G4double shieldingY = 5.*m;
  G4double shieldingZ = 10.*cm;
  G4double translation = -4.6* m;
  shielding = new G4Box("shielding",shieldingX/2.,shieldingY/2.,shieldingZ/2.);

  shieldingLog = new G4LogicalVolume(shielding,water,"shieldingLog",0,0,0);
  
  shieldingPhys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,translation + shieldingZ/2.),
            "shieldingPhys", shieldingLog, mother,false,0); 

  //Visualisation attributes
  G4Colour  red      (1.0,0.0,0.0);
  G4VisAttributes* shieldingVisAtt = new G4VisAttributes(red);
  shieldingVisAtt -> SetVisibility(true);
  shieldingVisAtt -> SetForceSolid(true);
  shieldingLog -> SetVisAttributes(shieldingVisAtt); 
}
