#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "RemSimVehicle1.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

RemSimVehicle1::RemSimVehicle1():layer1Phys(0)
{
 pMaterial = new RemSimMaterial();
}
RemSimVehicle1::~RemSimVehicle1()
{
  delete pMaterial;
}
void RemSimVehicle1::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
  
  G4double sizeX = 5.*m;
  G4double sizeY = 5.*m;

  G4Material* water = pMaterial ->GetMaterial("Water");  
  
  //layer of water
  G4double thick = 15. *cm;  

  G4Box* layer1 = new G4Box("layer1",sizeX/2.,sizeY/2.,thick/2.);

  G4LogicalVolume* layer1Log = new G4LogicalVolume(layer1,
                                                    water,
                                                    "layer1Log",
                                                    0,0,0);
  
  layer1Phys = new G4PVPlacement(0,
                                 G4ThreeVector(0.,0.,thick/2.),
                                 "layer1Phys", 
                                 layer1Log,
                                 motherVolume,
                                 false,
                                 0);

  G4Colour  blue   (0.,0.,1.);

  G4VisAttributes* layer1VisAtt = new G4VisAttributes(blue);
  layer1VisAtt -> SetVisibility(true);
  layer1VisAtt -> SetForceSolid(true);
  layer1Log -> SetVisAttributes(layer1VisAtt); 
}
void RemSimVehicle1::DestroyComponent()
{}
G4double  RemSimVehicle1::GetDensity()
{
  return 0;
}
void RemSimVehicle1::ChangeMaterial(G4String)
{}



