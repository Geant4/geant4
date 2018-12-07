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
/// \file electromagnetic/TestEm12/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorMaterial(nullptr),
   fAbsor(nullptr)
{
  // default parameter values
  fAbsorRadius = 3*cm;
  fNbOfLayers = 1;
  
  DefineMaterials();
  SetMaterial("G4_WATER");

  // create commands for interactive definition of the detector  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  G4NistManager* man = G4NistManager::Instance();
  
  man->FindOrBuildMaterial("G4_Al");
  man->FindOrBuildMaterial("G4_Si");
  man->FindOrBuildMaterial("G4_Fe");
  man->FindOrBuildMaterial("G4_Ge");
  man->FindOrBuildMaterial("G4_Gd");
  man->FindOrBuildMaterial("G4_W");
  man->FindOrBuildMaterial("G4_Pb");
  
  man->FindOrBuildMaterial("G4_AIR");
  man->FindOrBuildMaterial("G4_WATER");
  man->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
  
  ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Absorber
  //
  G4Sphere* 
  sAbsor = new G4Sphere("Absorber",                           //name
                     0., fAbsorRadius, 0., twopi, 0., pi);    //size

  fSpheres.push_back(sAbsor);

  G4LogicalVolume*
  lAbsor = new G4LogicalVolume(sAbsor,                        //solid
                                     fAbsorMaterial,          //material
                                    "Absorber");              //name
  fLVolumes.push_back(lAbsor);
                                   
  fAbsor = new G4PVPlacement(0,                         //no rotation
                             G4ThreeVector(),           //at (0,0,0)
                             lAbsor,                    //logical volume
                            "Absorber",                 //name
                             0,                         //mother  volume
                             false,                     //no boolean operation
                             0);                        //copy number

  // Layers
  //
  fLayerThickness = fAbsorRadius/fNbOfLayers;
                        
  for (G4int i=1; i<=fNbOfLayers; i++) {                           
    G4Sphere*
    sLayer = new G4Sphere("Layer", (i-1)*fLayerThickness, i*fLayerThickness,
                          0., twopi, 0., pi);

    fSpheres.push_back(sLayer);
                 
    G4LogicalVolume* 
    lLayer = new G4LogicalVolume(sLayer,                //shape
                                 fAbsorMaterial,        //material
                                 "Layer");              //name
    fLVolumes.push_back(lLayer);
                                 
             new G4PVPlacement(0,                      //no rotation
                               G4ThreeVector(),        //at (0,0,0)
                               lLayer,                 //logical volume
                               "Layer",                //name
                               lAbsor,                 //mother  volume
                               false,                  //no boolean operation
                               i);                     //copy number
                                                   
   }                                              

  PrintParameters();
    
  //
  //always return the root volume
  //  
  return fAbsor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Absorber is a sphere of " 
         << G4BestUnit(fAbsorRadius,"Length") << " radius of "
         << fAbsorMaterial->GetName() << " divided in " << fNbOfLayers 
         << " slices of " << G4BestUnit(fLayerThickness,"Length")
         << "\n \n" << fAbsorMaterial << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetRadius(G4double value)
{
  // geometry was already constructed - scale radii of all spheres
  if(fAbsor) {
    G4double scale = value/fAbsorRadius;
    for (auto solid : fSpheres) { 
      if(scale > 1.0) {
        solid->SetOuterRadius(solid->GetOuterRadius()*scale); 
        solid->SetInnerRadius(solid->GetInnerRadius()*scale);
      } else { 
        solid->SetInnerRadius(solid->GetInnerRadius()*scale);
        solid->SetOuterRadius(solid->GetOuterRadius()*scale); 
      }
    }
  }
  fAbsorRadius = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     

  if (pttoMaterial && pttoMaterial != fAbsorMaterial) {
    fAbsorMaterial = pttoMaterial;
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();

    // geometry was already constructed - only change material
    if(fAbsor) {
      for (auto lv : fLVolumes) { lv->SetMaterial(fAbsorMaterial); }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfLayers(G4int value)
{
  fNbOfLayers = value;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
