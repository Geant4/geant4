//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
// 
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorParameters.hh"
#include "CalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using namespace DetectorParameters;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  
  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4double calorXYSize = NofXYLayers * LayerXYSize;
  G4double calorZSize = NofZLayers * LayerZSize;
  G4double worldXYSize = 1.2 * calorXYSize;
  G4double worldZSize  = 1.2 * calorZSize; 
  
  // Get materials
  G4Material* defaultMaterial = G4Material::GetMaterial("G4_AIR");
  G4Material* absorberMaterial = G4Material::GetMaterial("G4_POLYETHYLENE");
  
  if ( ! defaultMaterial || ! absorberMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  G4VSolid* worldS 
    = new G4Box("World",           // its name
                 worldXYSize/2, worldXYSize/2, worldZSize/2); // its size
                         
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Calorimeter
  //  
  G4VSolid* calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorXYSize/2, calorXYSize/2, calorZSize/2); // its size
                         
  G4LogicalVolume* calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
   
  //                                 
  // Layer Z
  //
  G4VSolid* zlayerS 
    = new G4Box("ZLayer",           // its name
                 calorXYSize/2, calorXYSize/2, LayerZSize/2); //its size
                         
  G4LogicalVolume* zlayerLV
    = new G4LogicalVolume(
                 zlayerS,           // its solid
                 absorberMaterial,  // its material
                 "ZLayer");         // its name

  new G4PVReplica(
                 "ZLayer",          // its name
                 zlayerLV,          // its logical volume
                 calorLV,           // its mother
                 kZAxis,            // axis of replication
                 NofZLayers,        // number of replica
                 LayerZSize);       // witdth of replica
  
  //                                 
  // Layer YZ
  //
  G4VSolid* yzLayerS 
    = new G4Box("YZLayer",           // its name
                 calorXYSize/2, LayerXYSize/2, LayerZSize/2); //its size
                         
  G4LogicalVolume* yzLayerLV
    = new G4LogicalVolume(
                 yzLayerS,           // its solid
                 absorberMaterial,   // its material
                 "YZLayer");         // its name

  new G4PVReplica(
                 "YZLayer",          // its name
                 yzLayerLV,          // its logical volume
                 zlayerLV,           // its mother
                 kYAxis,             // axis of replication
                 NofXYLayers,        // number of replica
                 LayerXYSize);       // witdth of replica
  
  //                                 
  // Layer XYZ
  //
  G4VSolid* xyzLayerS 
    = new G4Box("XYZLayer",          // its name
                 LayerXYSize/2, LayerXYSize/2, LayerZSize/2); //its size
                         
  G4LogicalVolume* xyzLayerLV
    = new G4LogicalVolume(
                 xyzLayerS,          // its solid
                 absorberMaterial,   // its material
                 "XYZLayer");        // its name

  new G4PVReplica(
                 "XYZLayer",         // its name
                 xyzLayerLV,         // its logical volume
                 yzLayerLV,          // its mother
                 kXAxis,             // axis of replication
                 NofXYLayers,        // number of replica
                 LayerXYSize);       // witdth of replica
  
  //
  // print parameters
  //
  G4cout << "\n------------------------------------------------------------"
         << "\n---> The calorimeter is " << NofZLayers << " layers of: [ "
         << LayerZSize/mm << "mm of " << absorberMaterial->GetName() << " ] " 
         << "\n------------------------------------------------------------\n";
  
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  calorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //
  CalorimeterSD* absoSD = new CalorimeterSD("AbsorberSD");
  SetSensitiveDetector("XYZLayer", absoSD);

  // 
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
