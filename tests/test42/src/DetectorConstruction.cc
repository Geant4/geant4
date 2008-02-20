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
// $Id: DetectorConstruction.cc,v 1.7 2008-02-20 09:57:34 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// DetectorConstruction
//
// Created: 31.01.2003 V.Ivanchenko
//
// Modified:
// 09.12.2007 Adoptation of test42 (V.Grichine)
//
////////////////////////////////////////////////////////////////////////
// 

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "Tst42DetectorMessenger.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  logicTarget = 0;
  logicCheck  = 0;
  logicWorld  = 0;

  detectorMessenger  = new DetectorMessenger(this);
  fTst42DetMessenger = new Tst42DetectorMessenger(this);

  radius = 10.*cm;

  targetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_PbWo4");
  worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  // Prepare sensitive detectors
  checkSD = new CheckVolumeSD("checkSD");
  (G4SDManager::GetSDMpointer())->AddNewDetector( checkSD );
  targetSD = new TargetSD("targetSD");
  (G4SDManager::GetSDMpointer())->AddNewDetector( targetSD );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete detectorMessenger;
  delete fTst42DetMessenger;
}

/////////////////////////////////////////////////////////////////////
//
//


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructPhotoDetector();
}

///////////////////////////////////////////////////////////////////////
//
//

G4VPhysicalVolume* DetectorConstruction::ConstructPhotoDetector()
{
  if(fSetUp == "hadr01")
  {
    return Hadr01Construct();
  }
  else if(fSetUp == "simpleCMSPWO")
  {
    return CMSPWOsimpleConstruct();
  }
  else
  {
    return Hadr01Construct();
  }
}

///////////////////////////////////////////////////////////////////////////
//
// The first Hadr01 geometry

G4VPhysicalVolume* DetectorConstruction::Hadr01Construct()
{
  // Cleanup old geometry

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Sizes
  G4double checkR  = radius + mm;
  G4double worldR  = radius + cm;
  G4double targetZ = HistoManager::GetPointer()->Length()*0.5; 
  G4double checkZ  = targetZ + mm;
  G4double worldZ  = targetZ + cm;

  G4int nSlices    = HistoManager::GetPointer()->NumberOfSlices();
  G4double sliceZ  = targetZ/G4double(nSlices);

  //
  // World
  //
  // G4Tubs* solidW = new G4Tubs("World",0.,worldR,worldZ,0.,twopi);
  G4Box* solidW = new G4Box("World",worldR,worldR, worldZ);
  logicWorld = new G4LogicalVolume( solidW,worldMaterial,"World");
  G4VPhysicalVolume* world = new G4PVPlacement(0,G4ThreeVector(),
                                       logicWorld,"World",0,false,0);
  //
  // Check volume
  //
  // G4Tubs* solidC = new G4Tubs("Check",0.,checkR,checkZ,0.,twopi);
  G4Box* solidC = new G4Box("Check",checkR,checkR,checkZ);
  logicCheck = new G4LogicalVolume( solidC,worldMaterial,"World");
  //  G4VPhysicalVolume* physC = 
  new G4PVPlacement(0,G4ThreeVector(),logicCheck,"World",logicWorld,false,0);
  logicCheck->SetSensitiveDetector(checkSD);

  //
  // Target volume
  //
  // G4Tubs* solidA = new G4Tubs("Target",0.,radius,sliceZ,0.,twopi);
  G4Box* solidA = new G4Box("Target",radius,radius,sliceZ);
  logicTarget = new G4LogicalVolume( solidA,targetMaterial,"Target");
  logicTarget->SetSensitiveDetector(targetSD);

  G4double z = sliceZ - targetZ;

  for(G4int i=0; i<nSlices; i++) {
    // physC = 
    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,z),logicTarget,"Target",logicCheck,false,i);
    z += 2.0*sliceZ;
  }
  G4cout << "### Target consist of " << nSlices
         << " of " << targetMaterial->GetName() 
         << " disks with R(mm)= " << radius/mm
         << "  Width(mm)= " << 2.0*sliceZ/mm
         << "  Total Length(mm)= " << 2.0*targetZ/mm
         <<  "  ###" << G4endl;

//
// ------------ Generate & Add Material Properties Table ------------
//
  const G4int nEntries = 14;

  G4double PhotonEnergy[nEntries] =
            { 1.7712*eV,  1.8368*eV,  1.90745*eV, 1.98375*eV, 2.0664*eV, 
              2.15625*eV, 2.25426*eV, 2.3616*eV,  2.47968*eV, 2.61019*eV, 
              2.75521*eV, 2.91728*eV, 3.09961*eV, 3.30625*eV              };
//
// PbWO4
//
  G4double RefractiveIndex1[nEntries] =
            { 2.17728, 2.18025, 2.18357, 2.18753, 2.19285, 
              2.19813, 2.20441, 2.21337, 2.22328, 2.23619, 
              2.25203, 2.27381, 2.30282, 2.34666            };

 G4double Absorption1[nEntries] =
           { 666*cm, 666*cm, 666*cm, 666*cm, 666*cm, 
             666*cm, 605.455*cm, 512.308*cm, 444*cm, 333*cm, 
             246.667*cm, 195.882*cm, 195.882*cm, 158.571*cm         };

  G4double ScintilFast[nEntries] =
            { 0, 0, 0, 0, 0, 
              9, 23, 46, 72, 102, 
              121, 117, 84, 37                };

  G4double ScintilSlow[nEntries] =
            { 0, 0, 0, 0, 0, 
              9, 23, 46, 72, 102, 
              121, 117, 84, 37               };

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries);
  myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries);
  myMPT1->AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast,     nEntries);
  myMPT1->AddProperty("SLOWCOMPONENT",PhotonEnergy, ScintilSlow,     nEntries);

  G4double phYield = 200./MeV;
  phYield /= HistoManager::GetPointer()->GetPhotonBias();
  G4double birksConst = 0.;

  myMPT1->AddConstProperty("SCINTILLATIONYIELD",phYield);   // 2./MeV);
  myMPT1->AddConstProperty("BIRKSCONSTANT",birksConst);   
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT",10*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",36*ns);
  myMPT1->AddConstProperty("YIELDRATIO",1.0);

  targetMaterial->SetMaterialPropertiesTable(myMPT1);







  // colors
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* regWcolor = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  logicCheck->SetVisAttributes(regWcolor);

  G4VisAttributes* regCcolor = new G4VisAttributes(G4Colour(0., 0.3, 0.7));
  logicTarget->SetVisAttributes(regCcolor);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return world;
}

///////////////////////////////////////////////////////////////////////////
//
// Simple set-up with PbWO4 crystal with dimensions close to CMS ECAL

G4VPhysicalVolume* DetectorConstruction::CMSPWOsimpleConstruct()
{
  // Cleanup old geometry

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Sizes
  G4double checkR  = radius + mm;
  G4double worldR  = radius + cm;
  G4double targetZ = HistoManager::GetPointer()->Length()*0.5; 
  G4double checkZ  = targetZ + mm;
  G4double worldZ  = targetZ + cm;

  G4int nSlices    = HistoManager::GetPointer()->NumberOfSlices();
  G4double sliceZ  = targetZ/G4double(nSlices);

  //
  // World
  //
  // G4Tubs* solidW = new G4Tubs("World",0.,worldR,worldZ,0.,twopi);
  G4Box* solidW = new G4Box("World",worldR,worldR, worldZ);
  logicWorld = new G4LogicalVolume( solidW,worldMaterial,"World");
  G4VPhysicalVolume* world = new G4PVPlacement(0,G4ThreeVector(),
                                       logicWorld,"World",0,false,0);
  //
  // Check volume
  //
  // G4Tubs* solidC = new G4Tubs("Check",0.,checkR,checkZ,0.,twopi);
  G4Box* solidC = new G4Box("Check",checkR,checkR,checkZ);
  logicCheck = new G4LogicalVolume( solidC,worldMaterial,"World");
  //  G4VPhysicalVolume* physC = 
  new G4PVPlacement(0,G4ThreeVector(),logicCheck,"World",logicWorld,false,0);
  logicCheck->SetSensitiveDetector(checkSD);

  //
  // Target volume
  //
  // G4Tubs* solidA = new G4Tubs("Target",0.,radius,sliceZ,0.,twopi);
  G4Box* solidA = new G4Box("Target",radius,radius,sliceZ);
  logicTarget = new G4LogicalVolume( solidA,targetMaterial,"Target");
  logicTarget->SetSensitiveDetector(targetSD);

  G4double z = sliceZ - targetZ;

  for(G4int i=0; i<nSlices; i++) {
    // physC = 
    new G4PVPlacement(0,G4ThreeVector(0.0,0.0,z),logicTarget,"Target",logicCheck,false,i);
    z += 2.0*sliceZ;
  }
  G4cout << "### Target consist of " << nSlices
         << " of " << targetMaterial->GetName() 
         << " disks with R(mm)= " << radius/mm
         << "  Width(mm)= " << 2.0*sliceZ/mm
         << "  Total Length(mm)= " << 2.0*targetZ/mm
         <<  "  ###" << G4endl;

//
// ------------ Generate & Add Material Properties Table ------------
//
  const G4int nEntries = 14;

  G4double PhotonEnergy[nEntries] =
            { 1.7712*eV,  1.8368*eV,  1.90745*eV, 1.98375*eV, 2.0664*eV, 
              2.15625*eV, 2.25426*eV, 2.3616*eV,  2.47968*eV, 2.61019*eV, 
              2.75521*eV, 2.91728*eV, 3.09961*eV, 3.30625*eV              };
//
// PbWO4
//
  G4double RefractiveIndex1[nEntries] =
            { 2.17728, 2.18025, 2.18357, 2.18753, 2.19285, 
              2.19813, 2.20441, 2.21337, 2.22328, 2.23619, 
              2.25203, 2.27381, 2.30282, 2.34666            };

 G4double Absorption1[nEntries] =
           { 666*cm, 666*cm, 666*cm, 666*cm, 666*cm, 
             666*cm, 605.455*cm, 512.308*cm, 444*cm, 333*cm, 
             246.667*cm, 195.882*cm, 195.882*cm, 158.571*cm         };

  G4double ScintilFast[nEntries] =
            { 0, 0, 0, 0, 0, 
              9, 23, 46, 72, 102, 
              121, 117, 84, 37                };

  G4double ScintilSlow[nEntries] =
            { 0, 0, 0, 0, 0, 
              9, 23, 46, 72, 102, 
              121, 117, 84, 37               };

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries);
  myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries);
  myMPT1->AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast,     nEntries);
  myMPT1->AddProperty("SLOWCOMPONENT",PhotonEnergy, ScintilSlow,     nEntries);

  G4double phYield = 200./MeV;
  phYield /= HistoManager::GetPointer()->GetPhotonBias();
  G4double birksConst = 0.;

  myMPT1->AddConstProperty("SCINTILLATIONYIELD",phYield);   // 2./MeV);
  myMPT1->AddConstProperty("BIRKSCONSTANT",birksConst);   
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 0.1*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",0.3*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.8);

  targetMaterial->SetMaterialPropertiesTable(myMPT1);







  // colors
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* regWcolor = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  logicCheck->SetVisAttributes(regWcolor);

  G4VisAttributes* regCcolor = new G4VisAttributes(G4Colour(0., 0.3, 0.7));
  logicTarget->SetVisAttributes(regCcolor);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material && material != targetMaterial) {
    HistoManager::GetPointer()->SetTargetMaterial(material);
    targetMaterial = material;
    if(logicTarget) logicTarget->SetMaterial(targetMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material && material != worldMaterial) {
    worldMaterial = material;
    if(logicWorld) logicWorld->SetMaterial(worldMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetRadius(G4double val)  
{
  if(val > 0.0) {
    radius = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
