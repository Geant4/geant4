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
/// \file hadronic/Hadr02/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
/////////////////////////////////////////////////////////////////////////
//
// DetectorConstruction
//
// Created: 31.01.2003 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
// 

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4RunManager.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "TargetSD.hh"
#include "G4SDManager.hh"
#include "HistoManager.hh"
#include "G4NistManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fRadius(10.*cm),
   fTargetMaterial(0),
   fWorldMaterial(0),
   fTargetSD(0),
   fLogicTarget(0),
   fLogicWorld (0),
   fDetectorMessenger(0)
  
{
  fDetectorMessenger = new DetectorMessenger(this);

  fTargetMaterial 
    = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  fWorldMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  HistoManager::GetPointer()->SetTargetMaterial(fTargetMaterial);

  // Prepare sensitive detectors
  fTargetSD = new TargetSD("targetSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(fTargetSD);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Sizes
  G4double checkR  = fRadius + mm;
  G4double worldR  = fRadius + cm;
  G4double targetZ = HistoManager::GetPointer()->Length()*0.5; 
  G4double checkZ  = targetZ + mm;
  G4double worldZ  = targetZ + cm;

  G4int nSlices    = HistoManager::GetPointer()->NumberOfSlices();
  G4double sliceZ  = targetZ/G4double(nSlices);
 
  //
  // World
  G4Tubs* solidW = new G4Tubs("World",0.,worldR,worldZ,0.,twopi);
  fLogicWorld = new G4LogicalVolume( solidW,fWorldMaterial,"World");
  G4VPhysicalVolume* world = new G4PVPlacement(0,G4ThreeVector(),
                                       fLogicWorld,"World",0,false,0);
  //
  // Check volume
  //
  G4Tubs* solidC = new G4Tubs("Check",0.,checkR,checkZ,0.,twopi);
  G4LogicalVolume* logicCheck 
    = new G4LogicalVolume( solidC,fWorldMaterial,"World");
  new G4PVPlacement(0,G4ThreeVector(),logicCheck,"World",fLogicWorld,false,0);

  //
  // Target volume
  //
  G4Tubs* solidA = new G4Tubs("Target",0.,fRadius,sliceZ,0.,twopi);
  fLogicTarget = new G4LogicalVolume( solidA,fTargetMaterial,"Target");
  fLogicTarget->SetSensitiveDetector(fTargetSD);

  G4double z = sliceZ - targetZ;

  for(G4int i=0; i<nSlices; i++) {
    // physC = 
    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,z),fLogicTarget,
                      "Target",logicCheck,false,i);
    z += 2.0*sliceZ;
  }
  G4cout << "### Target consist of " << nSlices
         << " of " << fTargetMaterial->GetName() 
         << " disks with R(mm)= " << fRadius/mm
         << "  Width(mm)= " << 2.0*sliceZ/mm
         << "  Total Length(mm)= " << 2.0*targetZ/mm
         <<  "  ###" << G4endl;

  // colors
  G4VisAttributes zero = G4VisAttributes::GetInvisible();
  fLogicWorld->SetVisAttributes(zero);

  G4VisAttributes regWcolor(G4Colour(0.3, 0.3, 0.3));
  logicCheck->SetVisAttributes(regWcolor);

  G4VisAttributes regCcolor(G4Colour(0., 0.3, 0.7));
  fLogicTarget->SetVisAttributes(regCcolor);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);
 
  if (material && material != fTargetMaterial) {
    HistoManager::GetPointer()->SetTargetMaterial(material);
    fTargetMaterial = material;
    if(fLogicTarget) { fLogicTarget->SetMaterial(fTargetMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material && material != fWorldMaterial) {
    fWorldMaterial = material;
    if(fLogicWorld) { fLogicWorld->SetMaterial(fWorldMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetRadius(G4double val)  
{
  if(val > 0.0) {
    fRadius = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
