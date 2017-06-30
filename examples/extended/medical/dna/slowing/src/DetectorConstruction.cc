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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file medical/dna/slowing/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4NistManager.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction() :
    G4VUserDetectorConstruction(), fWaterMaterial(0)
{
  // create commands for interactive definition of the detector  
  fDetectorMessenger = new DetectorMessenger(this);

  //default tracking cut  
  fTrackingCut = 7.4*eV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()

{
  DefineMaterials();
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::DefineMaterials()
{

  // Water is defined from NIST material database
  G4NistManager * man = G4NistManager::Instance();

  G4Material * H2O = man->FindOrBuildMaterial("G4_WATER");

  /*
   If one wishes to test other density value for water material,
   one should use instead:
   G4Material * H2O = man->BuildMaterialWithNewDensity("G4_WATER_MODIFIED",
   "G4_WATER",1.100*g/cm3);

   Note: any string for "G4_WATER_MODIFIED" parameter is accepted
   and "G4_WATER" parameter should not be changed
   Both materials are created and can be selected from dna.mac
   */
  fWaterMaterial = H2O;

  //G4cout << "-> Density of water material (g/cm3)="
  // << fWaterMaterial->GetDensity()/(g/cm/cm/cm) << G4endl;

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{

  // WORLD VOLUME

  G4double worldSizeX = 1 * m;
  G4double worldSizeY = worldSizeX;
  G4double worldSizeZ = worldSizeX;

  G4Box* solidWorld = new G4Box("World",               //its name
      worldSizeX / 2, worldSizeY / 2, worldSizeZ / 2); //its size

  fLogicWorld = new G4LogicalVolume(solidWorld, //its solid
      fWaterMaterial, //its material
      "World"); //its name

  G4PVPlacement* physiWorld = new G4PVPlacement(
      0, //no rotation
      G4ThreeVector(), //at (0,0,0)
      fLogicWorld,
      "World", //its name
      0, //its mother  volume
      false, //no boolean operation
      0); //copy number

  // Visualization attributes
  G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0)); //White
  worldVisAtt->SetVisibility(true);
  fLogicWorld->SetVisAttributes(worldVisAtt);

  G4VisAttributes* worldVisAtt1 = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
  worldVisAtt1->SetVisibility(true);

  // Tracking cut
  fLogicWorld->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX,
    fTrackingCut));    

  PrintParameters();

  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // Search the material by its name   
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial
     (materialChoice);

  if (pttoMaterial)
  {
    fWaterMaterial = pttoMaterial;
    G4LogicalVolume* logicWorld =
        G4LogicalVolumeStore::GetInstance()->GetVolume("World");
    logicWorld->SetMaterial(fWaterMaterial);
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTrackingCut(G4double value)
{
  fTrackingCut = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters() const
{
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The tracking cut is set to " 
         << G4BestUnit(fTrackingCut,"Energy") << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
