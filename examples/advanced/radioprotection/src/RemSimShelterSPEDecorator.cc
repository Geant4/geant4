#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "RemSimShelterSPEDecorator.hh"
#include "RemSimDecorator.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "RemSimSensitiveDetector.hh"
#include "RemSimROGeometry.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"

RemSimShelterSPEDecorator::RemSimShelterSPEDecorator(RemSimVGeometryComponent* comp)
  : RemSimDecorator(comp)
{
 shelterSPEX = 5.*m;
 shelterSPEY = 5.*m;
 shelterSPEZ = 75.*cm;  
 translation = - 2.* m;
 pMaterial = new RemSimMaterial();

 shelterSPEVisAtt = 0;
}
RemSimShelterSPEDecorator::~RemSimShelterSPEDecorator()
{
  delete pMaterial;
}
void RemSimShelterSPEDecorator::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
  pMaterial -> DefineMaterials();
  RemSimDecorator::ConstructComponent(motherVolume);
  ConstructShelterSPE(motherVolume);
}

void RemSimShelterSPEDecorator::DestroyComponent()
{
 delete shelterSPEVisAtt;
 shelterSPEVisAtt = 0;
 
 delete shelterSPEPhys; 
 shelterSPEPhys = 0;

 delete shelterSPELog;
 shelterSPELog = 0;

 delete shelterSPE;
 shelterSPE = 0;
}
void RemSimShelterSPEDecorator::ConstructShelterSPE(G4VPhysicalVolume* motherVolume)
{
  // Geometry definition
  pMaterial -> DefineMaterials();

  G4Material* water = pMaterial -> GetMaterial("Water");
  
  shelterSPE = new G4Box("shelterSPE",shelterSPEX/2.,shelterSPEY/2.,shelterSPEZ/2.);

  shelterSPELog = new G4LogicalVolume(shelterSPE, water,
                                     "shelterSPELog",0,0,0);
  
  shelterSPEPhys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,translation + shelterSPEZ/2.),
            "shelterSPEPhys", shelterSPELog, motherVolume,false,0); 

  //Visualisation attributes
  G4Colour  red      (1.0,0.0,0.0);
  shelterSPEVisAtt = new G4VisAttributes(red);
  shelterSPEVisAtt -> SetVisibility(true);
  shelterSPEVisAtt -> SetForceSolid(true);
  shelterSPELog -> SetVisAttributes(shelterSPEVisAtt);  
  PrintDetectorParameters();
}

void RemSimShelterSPEDecorator::ChangeThickness(G4double)
{
  ;
}

void RemSimShelterSPEDecorator::ChangeMaterial(G4String materialName)
{ 
  G4cout << " This option is not available " << G4endl; 
}
void RemSimShelterSPEDecorator::PrintDetectorParameters()
{
  G4cout << "-----------------------------------------------------------------------"
         << G4endl
         << "the shelterSPE is a box whose thickness is: " << G4endl
         << ((shelterSPE -> GetZHalfLength())*2.)/cm
         << " cm along the Z axis"
         << G4endl
         << "material of the shelterSPE: "
         << shelterSPELog -> GetMaterial() -> GetName() <<G4endl
         << G4endl;
}
