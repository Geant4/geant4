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
// $Id: DetectorConstruction.cc,v 1.1 2006-06-02 19:00:00 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
/////////////////////////////////////////////////////////////////////////
//
// IION: Simple Phantom
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4RunManager.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "PhantomSD.hh"
#include "CheckVolumeSD.hh"
#include "G4SDManager.hh"
#include "HistoManager.hh"
#include "G4NistManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  logicA = 0;
  logicB = 0;
  logicW = 0;
  detectorMessenger = new DetectorMessenger(this);
  DefineMaterials();
  checkSD = new CheckVolumeSD("checkSD");
  (G4SDManager::GetSDMpointer())->AddNewDetector( checkSD );
  phantomSD = new PhantomSD("phantomSD");
  (G4SDManager::GetSDMpointer())->AddNewDetector( phantomSD );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);
  worldMaterial = man->FindOrBuildMaterial("G4_AIR");
  gapMaterial   = man->FindOrBuildMaterial("G4_AIR");
  absMaterial   = man->FindOrBuildMaterial("G4_WATER");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4double fac = 1.2;
  G4double fac2= 1.01;

  HistoManager* hManager = HistoManager::GetPointer();

  G4double radius = hManager->AbsRadius();
  G4double width  = hManager->AbsWidth();
  G4double gap    = hManager->GapWidth();
  G4int nAbs      = hManager->NumberOfAbs();

  G4double checkZ = 0.5*(width+gap)*nAbs;
  G4double worldZ = checkZ*fac;
  //
  // World
  //
  G4Tubs* solidW = new G4Tubs("World",0.,radius*fac,worldZ,0.,twopi);
  logicW = new G4LogicalVolume( solidW,worldMaterial,"World");
  G4VPhysicalVolume* world = new G4PVPlacement(0,G4ThreeVector(),
                                       logicW,"World",0,false,0);
  //
  // Check volume
  //
  G4Tubs* solidC = new G4Tubs("Check",0.,radius*fac2,checkZ*fac2,0.,twopi);
  G4LogicalVolume* logicC = new G4LogicalVolume( solidC,worldMaterial,"World");
  G4VPhysicalVolume* physC = new G4PVPlacement(0,G4ThreeVector(),
                                       logicC,"World",logicW,false,0);
  logicC->SetSensitiveDetector(checkSD);

  //
  // All gaps
  //
  G4Tubs* solidB = new G4Tubs("Gaps",0.,radius,checkZ,0.,twopi);
  logicB = new G4LogicalVolume( solidB,gapMaterial,"Gaps");
  new G4PVPlacement(0,G4ThreeVector(),logicB,"World",logicC,false,0);
  //
  // Phantom volume
  //
  G4Tubs* solidA = new G4Tubs("Phantom",0.,radius,width*0.5,0.,twopi);
  logicA = new G4LogicalVolume( solidA,absMaterial,"Phantom");
  logicA->SetSensitiveDetector(phantomSD);


  G4double z = -checkZ + 0.5*width + gap;
  genPosZ = -checkZ*fac2 - 0.001*mm;

  for(G4int i=0; i<nAbs; i++) {
    physC = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,z),
                                       logicA,"Phantom",logicB,false,i);
    z += width + gap;
  }
  G4cout << "### Phantom consist of " << nAbs
         << " disks of R(mm)= " << radius/mm
         << " disks of Width(mm)= " << width/mm
         << " with gap(mm)= " << gap/mm
         << " of " << absMaterial->GetName() <<  "  ###" << G4endl;

  // color regions

  G4VisAttributes* regWcolor = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  logicW->SetVisAttributes(regWcolor);

  G4VisAttributes* regCcolor = new G4VisAttributes(G4Colour(0., 0.3, 0.7));
  logicC->SetVisAttributes(regCcolor);

  //G4VisAttributes* regAcolor = new G4VisAttributes(G4Colour(0., 0.7, 0.3));
  //logicA->SetVisAttributes(regAcolor);
  logicA->SetVisAttributes(G4VisAttributes::Invisible);
  logicB->SetVisAttributes(G4VisAttributes::Invisible);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material) {
    absMaterial = material;
    if(logicA) logicA->SetMaterial(absMaterial);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material) {
    gapMaterial = material;
    if(logicB) logicB->SetMaterial(gapMaterial);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material) {
    worldMaterial = material;
    if(logicW) logicW->SetMaterial(worldMaterial);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
