#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "RemSimVehicle1.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

RemSimVehicle1::RemSimVehicle1():
  layer(0),layer1Log(0), layer1Phys(0),layer2Log(0), layer2Phys(0),
  layer1VisAtt(0), layer2VisAtt(0)
{
  pMaterial = new RemSimMaterial();;
}
RemSimVehicle1::~RemSimVehicle1()
{
  delete pMaterial;
}
void RemSimVehicle1::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
   pMaterial = new RemSimMaterial();
  
   //G4Material* lead = pMaterial->GetMaterial("Lead");
  G4Material* al = pMaterial ->GetMaterial("Aluminium"); 

  G4double layerX = 1.0*m;
  G4double layerY = 1.0*m;
  G4double layerZ = 10.0*cm;

  layer = new G4Box("layer",layerX/2.,layerY/2.,layerZ/2.);

  layer1Log = new G4LogicalVolume(layer,al,"layer1Log",0,0,0);
  
  layer1Phys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,-5.5*cm),
            "layer1Phys", layer1Log,motherVolume,false,0);
  
  layer2Log = new G4LogicalVolume(layer,al,"layer2Log",0,0,0);
  
  layer2Phys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,5.5*cm),
            "layer2Phys", layer2Log,motherVolume,false,0);
 
  //Visual Attributes
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  lblue   (0.0, 0.0, .75);

  G4VisAttributes* layer1VisAtt = new G4VisAttributes(lblue);
  layer1VisAtt -> SetVisibility(true);
  layer1VisAtt ->SetForceSolid(true);
  layer1Log -> SetVisAttributes(layer1VisAtt); 
 
  G4VisAttributes* layer2VisAtt = new G4VisAttributes(magenta);
  layer2VisAtt -> SetVisibility(true);
  layer2VisAtt ->SetForceSolid(true);
  layer2Log -> SetVisAttributes(layer2VisAtt); 
  }

void RemSimVehicle1::DestroyComponent()
{
  delete layer1VisAtt;layer1VisAtt=0;
  delete layer2VisAtt;layer2VisAtt=0;
  delete layer1Phys;layer1Phys=0;
  delete layer2Phys;layer2Phys=0;
  delete layer1Log;layer1Log=0;
  delete layer2Log;layer2Log=0;
  delete layer; layer=0;
}
