#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "RemSimMoonHabitat.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "RemSimDecorator.hh"
#include "RemSimAstronautDecorator.hh"
RemSimMoonHabitat::RemSimMoonHabitat()
{
 pMaterial = new RemSimMaterial(); 
 shelterPhys = 0;
}
RemSimMoonHabitat::~RemSimMoonHabitat()
{
  delete pMaterial;
}
void RemSimMoonHabitat::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
     pMaterial -> DefineMaterials();

     G4Material* moonMaterial = pMaterial->GetMaterial("moon") ;
     G4Material* air = pMaterial->GetMaterial("Air") ;
     G4double sizeX = 30.*m;
     G4double sizeY = 30.*m;
     G4double sizeZ = 15.*m;
     G4Box* moon = new G4Box("moon",sizeX,sizeY,sizeZ);
     
     G4LogicalVolume* moonLogical = new G4LogicalVolume(moon, 
                                                        moonMaterial,
                                                        "moon",0,0,0);

     G4VPhysicalVolume* moonPhys = new G4PVPlacement(0,
                                                     G4ThreeVector(0.,0.,sizeZ),                                                     "moon",moonLogical,
                                                     motherVolume,
                                                     false,0);

  G4double thick = 4.5 *m;
  G4Box* shelter = new G4Box("shelter",4.5/2.*m,10./2.*m,thick/2.);
  G4LogicalVolume* shelterLog = new G4LogicalVolume(shelter,
                                                    air,
                                                    "shelter",
                                                    0,0,0);
 
  G4double translation = - sizeZ + thick/2. + 0.5 *m;
 
  shelterPhys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,translation),
                                  "shelter",shelterLog, 
                                  moonPhys,false,0);

//     G4double roofThick = 3.*m;
//      G4Trd* roof = new G4Trd("roof",2.25*m,7.*m,5.*m, 8.*m,roofThick/2.);

//      G4LogicalVolume* roofLog = new G4LogicalVolume(roof,
//                                                     moon,
//                                                     "roof",
//                                                      0,0,0);

//     G4VPhysicalVolume* roofPhys = new G4PVPlacement(0,
//                                     G4ThreeVector(0.,0.,thick/2.- roofThick/2.),
//                                     "roof",roofLog, 
//                                     spacePhys,false,0);


  
    G4Colour  red      (1.0,0.0,0.0);
    G4Colour  blue   (0.,0.,1.);

    G4VisAttributes* moonVisAtt = new G4VisAttributes(red);
    moonVisAtt -> SetVisibility(true);
    moonVisAtt -> SetForceWireframe(true);
    moonLogical -> SetVisAttributes(moonVisAtt); 
    //roofLog -> SetVisAttributes(moonVisAtt);   
 
   G4VisAttributes* shelterVisAtt = new G4VisAttributes(blue);
    shelterVisAtt -> SetVisibility(true);
    shelterVisAtt -> SetForceWireframe(true);
    shelterLog -> SetVisAttributes(shelterVisAtt);    
}

void RemSimMoonHabitat::DestroyComponent()
{
}
 
G4VPhysicalVolume* RemSimMoonHabitat::GetShelter()
{
  return shelterPhys;
}
