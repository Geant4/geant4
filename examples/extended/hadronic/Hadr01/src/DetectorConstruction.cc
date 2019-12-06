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
/// \file hadronic/Hadr01/src/DetectorConstruction.cc
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
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
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
#include "CheckVolumeSD.hh"
#include "G4SDManager.hh"
#include "HistoManager.hh"
#include "G4NistManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fTargetMaterial(nullptr),
   fWorldMaterial(nullptr),
   fSolidW(nullptr),
   fSolidA(nullptr),
   fSolidC(nullptr),
   fLogicTarget(nullptr),
   fLogicCheck(nullptr),
   fLogicWorld(nullptr),
   fPhysWorld(nullptr),
   fInitialized(false)
{
  fDetectorMessenger = new DetectorMessenger(this);

  fRadius = 10.*cm;

  G4NistManager* nist = G4NistManager::Instance();
  fTargetMaterial = nist->FindOrBuildMaterial("G4_Al");
  fWorldMaterial  = nist->FindOrBuildMaterial("G4_Galactic");

  //
  // define battery material using Bugzilla 2175 data
  //
  G4Element* elH  = nist->FindOrBuildElement(1);
  G4Element* elLi = nist->FindOrBuildElement(3);
  G4Element* elC  = nist->FindOrBuildElement(6);
  G4Element* elO  = nist->FindOrBuildElement(8);
  G4Element* elAl = nist->FindOrBuildElement(13);
  G4Element* elTi = nist->FindOrBuildElement(22);
  G4Element* elCo = nist->FindOrBuildElement(27);
  G4Element* elCu = nist->FindOrBuildElement(29);
  G4Material* bat = new G4Material("Battery",2.165*CLHEP::g/CLHEP::cm3,8);
  bat->AddElement(elC,  0.19518445618745);
  bat->AddElement(elAl, 0.398);
  bat->AddElement(elTi, 0.02);
  bat->AddElement(elCu, 0.084);
  bat->AddElement(elLi, 0.0170098229419813);
  bat->AddElement(elCo, 0.144570016541753);
  bat->AddElement(elO,  0.134206611504321);
  bat->AddElement(elH,  0.0070290928244947);
  bat->GetIonisation()->SetMeanExcitationEnergy(144.88*eV);

  ComputeGeomParameters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

void DetectorConstruction::ComputeGeomParameters()
{
  // Sizes
  fCheckR  = fRadius + CLHEP::mm;
  fWorldR  = fRadius + CLHEP::cm;
  fTargetZ = HistoManager::GetPointer()->Length()*0.5; 
  fCheckZ  = fTargetZ + CLHEP::mm;
  fWorldZ  = fTargetZ + CLHEP::cm;

  fSlices  = HistoManager::GetPointer()->NumberOfSlices();
  fSliceZ  = fTargetZ/G4double(fSlices);
  if(fPhysWorld) {
    fSolidW->SetOuterRadius(fWorldR);
    fSolidW->SetZHalfLength(fWorldZ);
    fSolidA->SetOuterRadius(fRadius);
    fSolidA->SetZHalfLength(fSliceZ);
    fSolidC->SetOuterRadius(fCheckR);
    fSolidC->SetZHalfLength(fCheckZ);
  }
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if(fPhysWorld) { return fPhysWorld; }
  ComputeGeomParameters();

  //
  // World
  //
  fSolidW = new G4Tubs("World",0.,fWorldR,fWorldZ,0.,twopi);
  fLogicWorld = new G4LogicalVolume(fSolidW, fWorldMaterial, "World");
  fPhysWorld = new G4PVPlacement(nullptr,G4ThreeVector(0.,0.,0.),
                                 fLogicWorld,"World",nullptr,false,0);
  //
  // Check volume
  //
  fSolidC = new G4Tubs("Check",0.,fCheckR,fCheckZ,0.,twopi);
  fLogicCheck = new G4LogicalVolume(fSolidC, fWorldMaterial, "Check"); 
  new G4PVPlacement(nullptr,G4ThreeVector(),fLogicCheck,"Check",
                    fLogicWorld,false,0);

  //
  // Target volume
  //
  fSolidA = new G4Tubs("Target",0.,fRadius,fSliceZ,0.,twopi);
  fLogicTarget = new G4LogicalVolume(fSolidA,fTargetMaterial,"Target");

  G4double z = fSliceZ - fTargetZ;

  for(G4int i=0; i<fSlices; ++i) {
    // physC = 
    new G4PVPlacement(nullptr,G4ThreeVector(0.0,0.0,z),fLogicTarget,"Target",
                      fLogicCheck,false,i);
    z += 2.0*fSliceZ;
  }
  G4cout << "### Target consist of " << fSlices
         << " disks of " << fTargetMaterial->GetName() 
         << " with R(mm)= " << fRadius/mm
         << "  Width(mm)= " << 2.0*fSliceZ/mm
         << "  Total Length(mm)= " << 2.0*fTargetZ/mm
         <<  "  ###" << G4endl;

  // colors
  G4VisAttributes zero = G4VisAttributes::GetInvisible();
  fLogicWorld->SetVisAttributes(zero);

  G4VisAttributes regWcolor(G4Colour(0.3, 0.3, 0.3));
  fLogicCheck->SetVisAttributes(regWcolor);

  G4VisAttributes regCcolor(G4Colour(0., 0.3, 0.7));
  fLogicTarget->SetVisAttributes(regCcolor);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return fPhysWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  if (!fInitialized) {
    // Prepare sensitive detectors
    CheckVolumeSD* fCheckSD = new CheckVolumeSD("checkSD");
    (G4SDManager::GetSDMpointer())->AddNewDetector( fCheckSD );
    fLogicCheck->SetSensitiveDetector(fCheckSD);

    TargetSD* fTargetSD = new TargetSD("targetSD");
    (G4SDManager::GetSDMpointer())->AddNewDetector( fTargetSD );
    fLogicTarget->SetSensitiveDetector(fTargetSD);
    fInitialized = true;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material && material != fTargetMaterial) {
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
    ComputeGeomParameters();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
