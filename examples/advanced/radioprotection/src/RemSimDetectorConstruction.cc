//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: RemSimDetectorConstruction.cc,v 1.6 2004-05-19 09:29:34 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RemSimDetectorConstruction.hh"
#include "RemSimMaterial.hh"
#include "RemSimDetectorMessenger.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "RemSimVGeometryComponent.hh"
#include "RemSimVehicle1.hh"
#include "G4VisAttributes.hh"
#include "RemSimSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "RemSimROGeometry.hh"
#include "RemSimDecorator.hh"
#include "RemSimShieldingDecorator.hh"
RemSimDetectorConstruction::RemSimDetectorConstruction()
  :  experimentalHall_log(0),experimentalHall_phys(0),
     phantomPhys(0),detectorPhys(0)
{
 pMaterial = new RemSimMaterial();
 pVehicle = new RemSimVehicle1();
 messenger = new RemSimDetectorMessenger(this);
 decoratorValue ="Nothing";
 decorator = 0;  
}

RemSimDetectorConstruction::~RemSimDetectorConstruction()
{ 
  delete decorator;
  delete messenger;
  delete pVehicle;
  delete pMaterial;
}

G4VPhysicalVolume* RemSimDetectorConstruction::Construct()
{ 
  pMaterial -> DefineMaterials();
  G4Material* vacuum = pMaterial -> GetMaterial("Galactic") ;
  G4Material* water = pMaterial -> GetMaterial("Water");
 
  G4double expHall_x = 30.0*m;
  G4double expHall_y = 30.0*m;
  G4double expHall_z = 30.0*m;
  
  G4Box* experimentalHall_box
    = new G4Box("world",expHall_x,expHall_y,expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,
                                             vacuum,"world",0,0,0);
  
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),"world",
                                      experimentalHall_log,0,
                                      false,0);
  ConstructVolume();

  // Phantom
  G4Box* phantom = new G4Box("phantom",2.5*m,2.5*m,25.*cm);
  G4LogicalVolume* phantomLog = new G4LogicalVolume(phantom,
                                                    water,
                                                    "phantom",
                                                    0,0,0);
 
  phantomPhys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,156.*cm),
                                  "phantom",phantomLog, 
                                  experimentalHall_phys,false,0);
  // Detector
  G4Box* detector = new G4Box("detector",2.7*m,2.7*m,2.5*cm);
  G4LogicalVolume* detectorLog = new G4LogicalVolume(detector,
                                                     vacuum,
                                                     "detector",
                                                     0,0,0);
 
  detectorPhys = new G4PVPlacement(0,
                                   G4ThreeVector(0.,0.,183.5*cm),
                                   "detector",detectorLog, 
                                   experimentalHall_phys,false,0);

  G4Colour  lblue   (0.0, 0.0, .75); 
  G4VisAttributes* phantomVisAtt = new G4VisAttributes(lblue);
  phantomVisAtt -> SetVisibility(true);
  phantomVisAtt -> SetForceSolid(true);
  phantomLog -> SetVisAttributes(phantomVisAtt); 
 
  detectorLog -> SetVisAttributes (G4VisAttributes::Invisible);
 
 //Sensitive Detector  
 G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String sensitiveDetectorName = "AstronautSD";
  RemSimSensitiveDetector* sensitiveDetector = new  
                                 RemSimSensitiveDetector(sensitiveDetectorName);
 
  G4double phantomX= 5. *m;
  G4double phantomY= 5. *m;
  G4double phantomZ= 30.*cm;
  G4int VoxelNbAlongZ= 30;

  RemSimROGeometry* ROGeometry = new RemSimROGeometry(phantomX,
                                                      phantomY,
                                                      phantomZ,
						      VoxelNbAlongZ);
  ROGeometry->BuildROGeometry();
  sensitiveDetector->SetROgeometry(ROGeometry);
  SDman->AddNewDetector(sensitiveDetector);
  phantomLog ->SetSensitiveDetector(sensitiveDetector);
  phantomLog->SetUserLimits(new G4UserLimits(0.1*cm));

  return experimentalHall_phys;
}
void RemSimDetectorConstruction::ConstructVolume()
{
  pVehicle -> ConstructComponent(experimentalHall_phys);
}
void RemSimDetectorConstruction::AddShielding(G4String value)
{ 
  decoratorValue = value;
 
  if (decoratorValue == "On")
    { 
     decorator = new RemSimShieldingDecorator(pVehicle);
     decorator -> ConstructComponent(experimentalHall_phys); 
     G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
    } 

if (decoratorValue == "Off")
    {
      decorator = retrieveDecorator();
      decorator -> DestroyComponent(); 
      G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
    } 
}

RemSimDecorator* RemSimDetectorConstruction::retrieveDecorator()
{
  return decorator;
}
void RemSimDetectorConstruction::ChangeShieldingMaterial(G4String material)
{
 if (decoratorValue == "On")
    { 
     decorator = retrieveDecorator();
     decorator -> ChangeMaterial(material);
     G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
    }
 else G4cout<<" Define the shielding before!"<< G4endl;

}
void RemSimDetectorConstruction::ChangeShieldingThickness(G4double thick)
{
 if (decoratorValue == "On")
    { 
     decorator = retrieveDecorator();
     decorator -> ChangeThickness(thick);
     G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
    }
 else G4cout<<" Define the shielding before!"<< G4endl;
}
