#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "RemSimVehicle2.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"

RemSimVehicle2::RemSimVehicle2():
  layer1(0), layer2(0), layer3(0), layer4(0), layer5(0),
  layer1Log(0), layer2Log(0),layer3Log(0),layer4Log(0),layer5Log(0),  
  layer1Phys(0), layer2Phys(0), layer3Phys(0), layer4Phys(0), layer5Phys(0)
{
  pMaterial = new RemSimMaterial();
}
RemSimVehicle2::~RemSimVehicle2()
{
  delete pMaterial;
}
void RemSimVehicle2::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
   pMaterial = new RemSimMaterial();
  
  G4Material* nylon = pMaterial->GetMaterial("nylon");
  G4Material* mylar = pMaterial ->GetMaterial("mylar"); 
  G4Material* beta = pMaterial ->GetMaterial("beta"); 
  G4Material* kevlar = pMaterial ->GetMaterial("kevlar"); 
 
  G4double rmin = 3.60*cm;
  G4double rmax = 3.75 *m;
  G4double height = 7.5 *m;
  layer1 = new G4Tubs("layer1",rmin,rmax,height/2.,0.,360.);

  layer1Log = new G4LogicalVolume(layer1,beta,"betalayerLog",0,0,0);
  
  layer1Phys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,0.),
            "betalayer1Phys", layer1Log,motherVolume,false,0);
  
  rmin = 3.60*cm;
  rmax = 3.70 *m;
  height = 7.0 *m;
  layer2 = new G4Tubs("layer2",rmin,rmax,height/2.,0.,360.);

  layer2Log = new G4LogicalVolume(layer2,mylar,"mylarlayer2Log",0,0,0);
  
  layer2Phys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,0.),
            "mylarlayer2Phys", layer2Log,layer1Phys,false,0);
  
  rmin = 3.60*cm;
  rmax = 3.65 *m;
  height = 6.5 *m;
  layer3 = new G4Tubs("layer3",rmin,rmax,height/2.,0.,360.);

  layer3Log = new G4LogicalVolume(layer3,kevlar,"kevlarlayer3Log",0,0,0);
  
  layer3Phys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,0.),
            "kevlarlayer3Phys", layer3Log,layer2Phys,false,0);
  
  rmin = 3.60*cm;
  rmax = 3.60 *m;
  height = 6.0 *m;
  layer4 = new G4Tubs("layer4",rmin,rmax,height/2.,0.,360.);

  layer4Log = new G4LogicalVolume(layer4,nylon,"nylonlayer4Log",0,0,0);
  
  layer4Phys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,0.),
            "nylonlayer4Phys", layer4Log,layer3Phys,false,0);
  
  rmin = 3.60*cm;
  rmax = 3.65 *m;
  height = 5.8 *m;
  layer5 = new G4Tubs("layer5",rmin,rmax,height/2.,0.,360.);

  layer5Log = new G4LogicalVolume(layer5,nylon,"nylonlayer5Log",0,0,0);
  
  layer5Phys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,0.),
            "nylonlayer5Phys", layer5Log,layer4Phys,false,0);
 
  //Visual Attributes
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  lblue   (0.0, 0.0, .75);
  G4Colour  red      (1.0,0.0,0.0);
  G4Colour  green    (0.0,1.0,0.0);
  G4Colour  pink  (0.5,0.0,0.0);

  layer1VisAtt = new G4VisAttributes(lblue);
  layer1VisAtt -> SetVisibility(true);
  layer1VisAtt ->SetForceWireframe(true);
  layer1Log -> SetVisAttributes(layer1VisAtt); 
 
  layer2VisAtt = new G4VisAttributes(magenta);
  layer2VisAtt -> SetVisibility(true);
  layer2VisAtt ->SetForceWireframe(true);
  layer2Log -> SetVisAttributes(layer2VisAtt); 
 
  layer3VisAtt = new G4VisAttributes(red);
  layer3VisAtt -> SetVisibility(true);
  layer3VisAtt ->SetForceWireframe(true);
  layer3Log -> SetVisAttributes(layer3VisAtt); 
 
  layer4VisAtt = new G4VisAttributes(green);
  layer4VisAtt -> SetVisibility(true);
  layer4VisAtt ->SetForceWireframe(true);
  layer4Log -> SetVisAttributes(layer4VisAtt); 

  layer5VisAtt = new G4VisAttributes(pink);
  layer5VisAtt -> SetVisibility(true);
  layer5VisAtt ->SetForceWireframe(true);
  layer5Log -> SetVisAttributes(layer5VisAtt); 
 }

void RemSimVehicle2::DestroyComponent()
{
  delete layer1VisAtt;layer1VisAtt=0;
  delete layer2VisAtt;layer2VisAtt=0;
  delete layer3VisAtt;layer3VisAtt=0;
  delete layer4VisAtt;layer4VisAtt=0;
  delete layer5VisAtt;layer5VisAtt=0; 

  delete layer5Phys;layer5Phys=0;
  delete layer5Log;layer5Log=0;
  delete layer5;layer5=0;
  

  delete layer4Phys;layer4Phys=0;
  delete layer4Log;layer4Log=0;
  delete layer4;layer4=0;

  delete layer3Phys;layer3Phys =0;
  delete layer3Log;layer3Log =0;
  delete layer3;layer3 =0;

  delete layer2Phys;layer2Phys =0;
  delete layer2Log;layer2Log =0;
  delete layer2;layer2 =0;

  delete layer1Phys;layer1Phys =0;
  delete layer1Log;layer1Log =0;
  delete layer1;layer1 =0;
}

