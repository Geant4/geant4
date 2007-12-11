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
// $Id: DetectorConstruction.cc,v 1.3 2007-12-11 14:40:50 grichine Exp $
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

  detectorMessenger = new DetectorMessenger(this);

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
}

G4VPhysicalVolume* DetectorConstruction::Construct()
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
  const G4int nEntries = 32;

  G4double PhotonEnergy[nEntries] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
//
// PbWO4
//
  G4double RefractiveIndex1[nEntries] =
            { 2.3435, 2.344,  2.3445, 2.345,  2.3455,
              2.346,  2.3465, 2.347,  2.3475, 2.348,
              2.3485, 2.3492, 2.35,   2.3505, 2.351,
              2.3518, 2.3522, 2.3530, 2.3535, 2.354,
              2.3545, 2.355,  2.3555, 2.356,  2.3568,
              2.3572, 2.358,  2.3585, 2.359,  2.3595,
              2.36,   2.3608};

 G4double Absorption1[nEntries] =
           {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
           15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
           45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
           52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
           30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
           17.500*m, 14.500*m };

  G4double ScintilFast[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };
  G4double ScintilSlow[nEntries] =
            { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
              7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
              3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
              4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
              7.00, 6.00, 5.00, 4.00 };

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries);
  myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries);
  myMPT1->AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast,     nEntries);
  myMPT1->AddProperty("SLOWCOMPONENT",PhotonEnergy, ScintilSlow,     nEntries);

  myMPT1->AddConstProperty("SCINTILLATIONYIELD",2./MeV);
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
