#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "RemSimAstronautDecorator.hh"
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
#include "G4UserLimits.hh"

RemSimAstronautDecorator::RemSimAstronautDecorator(RemSimVGeometryComponent* comp)
  : RemSimDecorator(comp),phantom(0), phantomLog(0), phantomPhys(0)
{ 
 phantomZ = 30.*cm;
 motherAstronaut = 0;
 flag = false;
}
RemSimAstronautDecorator::~RemSimAstronautDecorator()
{
}
void RemSimAstronautDecorator::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
  RemSimDecorator::ConstructComponent(motherVolume);
  if (flag == false) ConstructAstronaut(motherVolume);
  else ConstructAstronaut(motherAstronaut);
}

void RemSimAstronautDecorator::DestroyComponent()
{

}
void RemSimAstronautDecorator::ConstructAstronaut(G4VPhysicalVolume* motherVolume)
{
  // Astronaut definition: Box of water
  // The Astronaut is the sensitive detector
  RemSimMaterial*  pMaterial = new RemSimMaterial();

  G4Material* water = pMaterial -> GetMaterial("Water");
  delete pMaterial;

  G4double phantomX = 5. *m;
  G4double phantomY = 5. *m;
 
  phantom = new G4Box("phantom",phantomX/2.,phantomY/2.,phantomZ/2.);

  phantomLog = new G4LogicalVolume(phantom,
                                   water,
                                   "phantom",
                                   0,0,0);
 
  phantomPhys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,0.),
                                  "phantom",phantomLog, 
                                  motherVolume,false,0);

  G4Colour  lblue   (0.0, 0.0,.75); 
  G4VisAttributes* phantomVisAtt = new G4VisAttributes(lblue);
  phantomVisAtt -> SetVisibility(true);
  phantomVisAtt -> SetForceSolid(true);
  phantomLog -> SetVisAttributes(phantomVisAtt); 
 
 //Sensitive Detector  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String sensitiveDetectorName = "AstronautSD";
  
  RemSimSensitiveDetector* sensitiveDetector = new  
                                 RemSimSensitiveDetector(sensitiveDetectorName);
  G4int VoxelNbAlongZ = 30;
 
  G4double translation = 0;
    
  if (flag == true) 
    {
      G4double thickShelter = 4.5 *m;
      translation = 0.5*m + thickShelter/2.;
   }
 
  RemSimROGeometry* ROGeometry = new RemSimROGeometry(phantomX,
                                                      phantomY,
                                                      phantomZ,
						      VoxelNbAlongZ, 
                                                      translation);
  ROGeometry -> BuildROGeometry();
  sensitiveDetector -> SetROgeometry(ROGeometry);
  SDman->AddNewDetector(sensitiveDetector);
  phantomLog -> SetSensitiveDetector(sensitiveDetector);
  phantomLog -> SetUserLimits(new G4UserLimits(0.1*cm));
}

void RemSimAstronautDecorator::ChangeThickness(G4double thick)
{
  G4cout << "It is not possible to change the sizes of the Astronaut"
         << G4endl;
}

void RemSimAstronautDecorator::ChangeMaterial(G4String materialName)
{ 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);     
   if (pttoMaterial)
    {
      phantomLog -> SetMaterial(pttoMaterial);
      PrintDetectorParameters();
    } 
  else
    G4cout << "WARNING: material '" << materialName
           << "' not available!" << G4endl;         
}
void RemSimAstronautDecorator::PrintDetectorParameters()
{
  G4cout << "-----------------------------------------------------------------------"
         << G4endl
         << "the astronaut is a box whose thickness is: " << G4endl
         << phantomZ/cm
         << " cm along the Z axis"
         << G4endl
         << "material of the astronaut: "
         << phantomLog -> GetMaterial() -> GetName() <<G4endl
         << G4endl;
}
void RemSimAstronautDecorator::ChangeMother(G4VPhysicalVolume* mother)
{
 motherAstronaut = mother; 
 flag = true;                      
}
