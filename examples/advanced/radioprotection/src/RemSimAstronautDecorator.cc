#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
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

RemSimAstronautDecorator::
RemSimAstronautDecorator(RemSimVGeometryComponent* comp)
  : RemSimDecorator(comp)
{
  
}
RemSimAstronautDecorator::~RemSimAstronautDecorator()
{
  ;
}
void RemSimAstronautDecorator::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
   RemSimDecorator::ConstructComponent(motherVolume);
  ConstructAstronaut(motherVolume);
}

void RemSimAstronautDecorator::DestroyComponent()
{
  ;
}
void RemSimAstronautDecorator::ConstructAstronaut(G4VPhysicalVolume* mother)
{
  // Materials for the astronaut composition
  RemSimMaterial* pMaterial = new RemSimMaterial();
  G4Material* water = pMaterial -> GetMaterial("Water");
  delete pMaterial;

  // Geometry definition
  G4double astronautX = 10.*cm;
  G4double astronautY = 10* cm;
  G4double astronautZ = 10* cm;

  astronaut = new G4Box("astronaut",astronautX/2.,astronautY/2.,astronautZ/2.);

  astronautLog = new G4LogicalVolume(astronaut,water,"astronautLog",0,0,0);
  
  astronautPhys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,40.*cm),
            "astronautPhys", astronautLog,mother,false,0); 

  //Visualisation attributes
  G4Colour color(0.5,0.5,0.);
  G4VisAttributes* astronautVisAtt = new G4VisAttributes(color);
  astronautVisAtt -> SetVisibility(true);
  astronautVisAtt -> SetForceSolid(true);
  astronautLog -> SetVisAttributes(astronautVisAtt); 

  // sensitivity
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String sensitiveDetectorName = "AstronautSD";
  RemSimSensitiveDetector* sensitiveDetector = new  
                                 RemSimSensitiveDetector(sensitiveDetectorName);
  G4int VoxelNbAlongX = 100;
  G4int VoxelNbAlongY = 100;
  G4int VoxelNbAlongZ = 100;

  RemSimROGeometry* ROGeometry = new RemSimROGeometry(astronautX,
                                                            astronautY,
                                                            astronautZ,
							    VoxelNbAlongX,
							    VoxelNbAlongY, 
							    VoxelNbAlongZ);
  ROGeometry->BuildROGeometry();
  sensitiveDetector->SetROgeometry(ROGeometry);

  SDman->AddNewDetector(sensitiveDetector);
  astronautLog ->SetSensitiveDetector(sensitiveDetector);
}
