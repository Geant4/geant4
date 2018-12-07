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
/// \file electromagnetic/TestEm9/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//
/////////////////////////////////////////////////////////////////////////
//
// TestEm9: Crystal calorimeter
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

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4TransportationManager.hh"

#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4NistManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fCalMaterial(0),
    fVertMaterial(0),
    fAbsMaterial(0),
    fWorldMaterial(0),
    fYorkMaterial(0),
    fLogicWorld(0),
    fLogicCal(0),
    fLogicA1(0),
    fLogicA2(0),
    fLogicA3(0),
    fLogicA4(0),
    fVertexRegion(0),
    fMuonRegion(0),
    fVertexDetectorCuts(0),
    fMuonDetectorCuts(0),
    fDetectorMessenger(0)
{
  fDetectorMessenger = new DetectorMessenger(this);

  fEcalLength   = 36.*cm;
  fEcalWidth    = 6.*cm;
  fVertexLength = 3.*cm;
  fPadLength    = 0.1*mm;
  fPadWidth     = 0.02*mm;
  fAbsLength    = 2.*mm;
  fWorldZ       = 0.0;
  fLogicWorld   = 0;
  fLogicCal     = 0;
  fLogicA1      = 0;
  fLogicA2      = 0;
  fLogicA3      = 0;
  fLogicA4      = 0;
  fVertexRegion = 0;
  fMuonRegion   = 0;

  DefineMaterials();
  fVertexDetectorCuts = new G4ProductionCuts();
  fMuonDetectorCuts   = new G4ProductionCuts();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
  delete fVertexDetectorCuts;
  delete fMuonDetectorCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Default materials

  G4NistManager* man = G4NistManager::Instance();
  //  man->SetVerbose(1);
  fWorldMaterial = man->FindOrBuildMaterial("G4_AIR");
  fAbsMaterial   = man->FindOrBuildMaterial("G4_Al");
  fVertMaterial  = man->FindOrBuildMaterial("G4_Si");
  fYorkMaterial  = man->FindOrBuildMaterial("G4_Fe");
  fCalMaterial   = man->FindOrBuildMaterial("G4_CESIUM_IODIDE");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry

  G4GeometryManager::GetInstance()->OpenGeometry();

  if(G4NistManager::Instance()->GetVerbose() > 0)
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  if(fVertexRegion) {
    delete fVertexRegion;
    delete fMuonRegion;
  }
  fVertexRegion = new G4Region("VertexDetector");
  fVertexRegion->SetProductionCuts(fVertexDetectorCuts);

  fMuonRegion   = new G4Region("MuonDetector");
  fMuonRegion->SetProductionCuts(fMuonDetectorCuts);

  G4SolidStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4PhysicalVolumeStore::GetInstance()->Clean();

  if(fVertexLength < fPadLength*5.0) fVertexLength = fPadLength*5.0;
  G4double gap    = 0.01*mm;
  G4double biggap = 2.*cm;
  G4double york   = 10.*cm;

  fWorldZ = 2.*fVertexLength + 3.*fAbsLength + 0.5*(fEcalLength + york) + biggap*2.;

  G4double worldX = fEcalWidth*3.0;
  G4double vertexZ= -fWorldZ + fVertexLength*2.0 + fAbsLength     + biggap;
  G4double absZ2  = -fWorldZ + fVertexLength*4.0 + fAbsLength*3.5 + biggap;
  G4double ecalZ  = -fWorldZ + fVertexLength*4.0 + fAbsLength*4.0 + fEcalLength*0.5 
    + 2.*biggap;
  G4double yorkZ  = -fWorldZ + fVertexLength*4.0 + fAbsLength*5.0 + fEcalLength
    + york*0.5 + 3.*biggap;

  //
  // World
  //
  G4Box* solidW = new G4Box("World",worldX,worldX,fWorldZ);
  fLogicWorld = new G4LogicalVolume( solidW,fWorldMaterial,"World");
  G4VPhysicalVolume* world = new G4PVPlacement(0,G4ThreeVector(),
                                               "World",fLogicWorld,0,false,0);

  //
  // Ecal
  //
  G4Box* solidE = new G4Box("VolE",worldX,worldX,fEcalLength*0.5 + gap);
  G4LogicalVolume* logicECal = 
    new G4LogicalVolume( solidE,fWorldMaterial,"VolE");
  G4VPhysicalVolume* physE = new G4PVPlacement(0,G4ThreeVector(0.,0.,ecalZ),
                                               "VolE",logicECal,world,false,0);

  G4Box* solidC = new G4Box("Ecal",fEcalWidth*0.5,fEcalWidth*0.5,fEcalLength*0.5);
  fLogicCal = new G4LogicalVolume( solidC,fCalMaterial,"Ecal");

  G4cout << "Ecal is " << G4BestUnit(fEcalLength,"Length")
         << " of " << fCalMaterial->GetName() << G4endl;

  // Crystals

  G4double x0 = -(fEcalWidth + gap)*2.0;
  G4double y  = x0;
  G4double x;
  G4int k = 0;
  G4int i,j;

  for (i=0; i<5; i++) {
    x  = x0;
    for (j=0; j<5; j++) {

      new G4PVPlacement(0,G4ThreeVector(x,y,0.),"Ecal",fLogicCal,
                                    physE,false,k);
      k++;
      x += fEcalWidth + gap;
    }
    y += fEcalWidth + gap;
  }

  //Absorber

  G4Box* solidA = new G4Box("Abso",worldX,worldX,fAbsLength*0.5);
  fLogicA2 = new G4LogicalVolume( solidA,fAbsMaterial,"Abs2");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,absZ2),
                         "Abs2",fLogicA2,world,false,0);

  G4cout << "Absorber is " << G4BestUnit(fAbsLength,"Length")
         << " of " << fAbsMaterial->GetName() << G4endl;

  //York

  G4Box* solidYV = new G4Box("VolY",worldX,worldX,york*0.5+fAbsLength);
  G4LogicalVolume* logicYV = 
    new G4LogicalVolume( solidYV,fYorkMaterial,"VolY");
  G4VPhysicalVolume* physYV = new G4PVPlacement(0,G4ThreeVector(0.,0.,yorkZ),
                                                "VolY",logicYV,world,false,0);

  G4Box* solidY = new G4Box("York",worldX,worldX,york*0.5);
  G4LogicalVolume* logicY = 
    new G4LogicalVolume( solidY,fYorkMaterial,"York");
  new G4PVPlacement(0,G4ThreeVector(),
                         "York",logicY,physYV,false,0);

  fLogicA3 = new G4LogicalVolume( solidA,fAbsMaterial,"Abs3");
  fLogicA4 = new G4LogicalVolume( solidA,fAbsMaterial,"Abs4");

  new G4PVPlacement(0,G4ThreeVector(0.,0.,-(york+fAbsLength)*0.5),
                         "Abs3",fLogicA3,physYV,false,0);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,(york+fAbsLength)*0.5),
                         "Abs4",fLogicA4,physYV,false,0);

  //Vertex volume
  G4Box* solidVV = new G4Box("VolV",worldX,worldX,fVertexLength*2.+fAbsLength+gap);
  G4LogicalVolume* logicVV = 
    new G4LogicalVolume( solidVV,fWorldMaterial,"VolV");
  G4VPhysicalVolume* physVV = new G4PVPlacement(0,G4ThreeVector(0.,0.,vertexZ),
                                                "VolV",logicVV,world,false,0);

  //Absorber
  fLogicA1 = new G4LogicalVolume( solidA,fAbsMaterial,"Abs1");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,fVertexLength*2.-fAbsLength*0.5),
                         "Abs1",fLogicA1,physVV,false,0);

  //Vertex
  G4double vertWidth = fEcalWidth/5.;
  G4int npads = (G4int)(vertWidth/fPadWidth);
  //G4cout << " vertWidth= " << vertWidth << " padWidth= " << padWidth 
  //         << " npads= " << npads << G4endl;
  // insure beam to hit a middle of central pad
  npads = (npads/2)*2 + 1;
  x0 = -0.5*(fPadWidth + vertWidth);
  G4double x1 = 0.5*vertWidth + gap; 
  G4double z  = -(fVertexLength+fAbsLength);

  G4Box* solidVD = new G4Box("VertDet",x1,fEcalWidth*0.5+gap,fPadLength*0.5);
  G4LogicalVolume* logicVD = 
    new G4LogicalVolume( solidVD,fVertMaterial,"VertDet");
  logicVD->SetSolid(solidVD);

  G4Box* solidV = new G4Box("Vert",fPadWidth*0.5,fEcalWidth*0.5,fPadLength*0.5);
  G4LogicalVolume* logicV = new G4LogicalVolume( solidV,fVertMaterial,"Vert");

  for (i=0; i<3; i++) {
    new G4PVPlacement(0,G4ThreeVector(0.,0.,z),"VertDet",logicVD,
                                    physVV,false,i);
    z += fVertexLength;
  }
  x = x0;

  for (j=0; j<npads; j++) {

    new G4PVPlacement(0,G4ThreeVector(x,0.,0.),logicV,"Vert",logicVD,
                                    false,k);
    x += fPadWidth;
  }

  G4cout << "Vertex is " << G4BestUnit(fVertexLength,"Length")
         << " of 3 layers of Si of " << G4BestUnit(fPadLength,"Length")
         << " npads= " << npads
         << G4endl;

  // Define region for the vertex detector
  fVertexRegion->AddRootLogicalVolume(logicVV);
  fVertexRegion->AddRootLogicalVolume(fLogicA3);

  // Define region for the muon detector
  fMuonRegion->AddRootLogicalVolume(logicYV);

  // color regions
  logicVV-> SetVisAttributes(G4VisAttributes::GetInvisible());
  logicV-> SetVisAttributes(G4VisAttributes::GetInvisible());
  logicECal-> SetVisAttributes(G4VisAttributes::GetInvisible());
  logicYV-> SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VisAttributes* regWcolor = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  fLogicWorld->SetVisAttributes(regWcolor);

  G4VisAttributes* regVcolor = new G4VisAttributes(G4Colour(0., 0.3, 0.7));
  logicVD->SetVisAttributes(regVcolor);

  G4VisAttributes* regCcolor = new G4VisAttributes(G4Colour(0., 0.7, 0.3));
  fLogicCal->SetVisAttributes(regCcolor);

  G4VisAttributes* regAcolor = new G4VisAttributes(G4Colour(1., 0.5, 0.5));
  fLogicA1->SetVisAttributes(regAcolor);
  fLogicA2->SetVisAttributes(regAcolor);
  fLogicA3->SetVisAttributes(regAcolor);
  fLogicA4->SetVisAttributes(regAcolor);

  G4VisAttributes* regMcolor = new G4VisAttributes(G4Colour(1., 1., 0.));
  logicY->SetVisAttributes(regMcolor);

  // always return world
  G4cout << "### New geometry is constructed" << G4endl;
 
  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetEcalMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if (pttoMaterial && pttoMaterial != fCalMaterial) {
    fCalMaterial = pttoMaterial;
    if(fLogicCal) {
      fLogicCal->SetMaterial(fCalMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if (pttoMaterial && pttoMaterial != fAbsMaterial) {
    fAbsMaterial = pttoMaterial;
    if(fLogicA1) {
      fLogicA1->SetMaterial(fAbsMaterial);
      fLogicA2->SetMaterial(fAbsMaterial);
      fLogicA3->SetMaterial(fAbsMaterial);
      fLogicA4->SetMaterial(fAbsMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetEcalLength (G4double val)   
{
  if(val > 0.0) {
    fEcalLength = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetEcalWidth  (G4double val)   
{
  if(val > 0.0) {
    fEcalWidth = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetVertexLength (G4double val) 
{
  if(val > 0.0) {
    fVertexLength = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPadLength  (G4double val)   
{
  if(val > 0.0) {
    fPadLength = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPadWidth  (G4double val)    
{
  if(val > 0.0) {
    fPadWidth = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsLength(G4double val)     
{
  if(val > 0.0) {
    fAbsLength = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
