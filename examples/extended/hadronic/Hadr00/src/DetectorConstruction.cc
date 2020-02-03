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
/// \file hadronic/Hadr00/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
<<<<<<< HEAD
// $Id: DetectorConstruction.cc 77254 2013-11-22 10:08:02Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
/////////////////////////////////////////////////////////////////////////
//
// DetectorConstruction
//
// Created: 20.06.2008 V.Ivanchenko
//
// Modified:
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

#include "G4NistManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fTargetMaterial(0), fWorldMaterial(0),
  fLogicTarget(0), fLogicWorld(0),
  fDetectorMessenger(0)
{
  fDetectorMessenger = new DetectorMessenger(this);

  fRadius = 5.*cm;
  fLength = 10.*cm;

  fTargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  fWorldMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
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
  G4double worldR  = fRadius + cm;
  G4double worldZ  = fLength + cm;

  //
  // World
  //
  G4Tubs* solidW = new G4Tubs("World", 0., worldR, 0.5*worldZ, 0., twopi);
  fLogicWorld = new G4LogicalVolume(solidW, fWorldMaterial,"World");
  G4VPhysicalVolume* world = new G4PVPlacement(0, G4ThreeVector(),
                                               fLogicWorld, "World", 
                                               0, false, 0);

  //
  // Target volume
  //
  G4Tubs* solidA = new G4Tubs("Target", 0., fRadius, 0.5*fLength, 0.,twopi);
  fLogicTarget = new G4LogicalVolume( solidA, fTargetMaterial, "Target");
  new G4PVPlacement(0, G4ThreeVector(), fLogicTarget, "Target",
                    fLogicWorld, false, 0);

  G4cout << "### Target consist of " 
         << fTargetMaterial->GetName() 
         << " disks with R(mm)= " << fRadius/mm
         << "  fLength(mm)= " << fLength/mm
         <<  "  ###" << G4endl;

  // colors
  fLogicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* regCcolor = new G4VisAttributes(G4Colour(0., 0.3, 0.7));
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
    fTargetMaterial = material;
    if(fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
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
    if(fLogicWorld) fLogicWorld->SetMaterial(fWorldMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

void DetectorConstruction::SetTargetRadius(G4double val)  
{
  if(val > 0.0 && val != fRadius) {
    fRadius = val;
      G4RunManager::GetRunManager()->ReinitializeGeometry();
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetLength(G4double val)  
{
  if(val > 0.0 && val != fLength) {
    fLength = val;
      G4RunManager::GetRunManager()->ReinitializeGeometry();
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
