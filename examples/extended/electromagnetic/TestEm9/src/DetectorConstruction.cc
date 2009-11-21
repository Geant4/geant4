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
//
// $Id: DetectorConstruction.cc,v 1.14 2009-11-21 17:28:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction()
{
  detectorMessenger = new DetectorMessenger(this);

  ecalLength   = 36.*cm;
  ecalWidth    = 6.*cm;
  vertexLength = 3.*cm;
  padLength    = 0.1*mm;
  padWidth     = 0.02*mm;
  absLength    = 2.*mm;
  vertexRegion = 0;
  muonRegion   = 0;
  logicWorld   = 0;
  logicCal     = 0;
  logicA1      = 0;
  DefineMaterials();
  vertexDetectorCuts = new G4ProductionCuts();
  muonDetectorCuts   = new G4ProductionCuts();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete detectorMessenger;
  delete vertexDetectorCuts;
  delete muonDetectorCuts;
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
  worldMaterial = man->FindOrBuildMaterial("G4_AIR");
  absMaterial   = man->FindOrBuildMaterial("G4_Al");
  vertMaterial  = man->FindOrBuildMaterial("G4_Si");
  yorkMaterial  = man->FindOrBuildMaterial("G4_Fe");
  calMaterial   = man->FindOrBuildMaterial("G4_CESIUM_IODIDE");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry

  G4GeometryManager::GetInstance()->OpenGeometry();

  if(G4NistManager::Instance()->GetVerbose() > 0)
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  if(vertexRegion) {
    delete vertexRegion;
    delete muonRegion;
  }
  vertexRegion = new G4Region("VertexDetector");
  vertexRegion->SetProductionCuts(vertexDetectorCuts);

  muonRegion   = new G4Region("MuonDetector");
  muonRegion->SetProductionCuts(muonDetectorCuts);

  G4SolidStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4PhysicalVolumeStore::GetInstance()->Clean();

  if(vertexLength < padLength*5.0) vertexLength = padLength*5.0;
  G4double gap    = 0.01*mm;
  G4double biggap = 2.*cm;
  G4double york   = 10.*cm;

  worldZ = 2.*vertexLength + 3.*absLength + 0.5*(ecalLength + york) + biggap*2.;

  G4double worldX = ecalWidth*3.0;
  G4double vertexZ= -worldZ + vertexLength*2.0 + absLength     + biggap;
  G4double absZ2  = -worldZ + vertexLength*4.0 + absLength*3.5 + biggap;
  G4double ecalZ  = -worldZ + vertexLength*4.0 + absLength*4.0 + ecalLength*0.5 
    + 2.*biggap;
  G4double yorkZ  = -worldZ + vertexLength*4.0 + absLength*5.0 + ecalLength
    + york*0.5 + 3.*biggap;

  //
  // World
  //
  G4Box* solidW = new G4Box("World",worldX,worldX,worldZ);
  logicWorld = new G4LogicalVolume( solidW,worldMaterial,"World");
  G4VPhysicalVolume* world = new G4PVPlacement(0,G4ThreeVector(),
					       "World",logicWorld,0,false,0);

  //
  // Ecal
  //
  G4Box* solidE = new G4Box("VolE",worldX,worldX,ecalLength*0.5 + gap);
  logicECal = new G4LogicalVolume( solidE,worldMaterial,"VolE");
  G4VPhysicalVolume* physE = new G4PVPlacement(0,G4ThreeVector(0.,0.,ecalZ),
					       "VolE",logicECal,world,false,0);

  G4Box* solidC = new G4Box("Ecal",ecalWidth*0.5,ecalWidth*0.5,ecalLength*0.5);
  logicCal = new G4LogicalVolume( solidC,calMaterial,"Ecal");

  G4cout << "Ecal is " << G4BestUnit(ecalLength,"Length")
	 << " of " << calMaterial->GetName() << G4endl;

  // Crystals

  G4double x0 = -(ecalWidth + gap)*2.0;
  G4double y  = x0;
  G4double x;
  G4int k = 0;
  G4VPhysicalVolume* pv;
  G4int i,j;

  for (i=0; i<5; i++) {
    x  = x0;
    for (j=0; j<5; j++) {

      pv = new G4PVPlacement(0,G4ThreeVector(x,y,0.),"Ecal",logicCal,
                                    physE,false,k);
      k++;
      x += ecalWidth + gap;
    }
    y += ecalWidth + gap;
  }

  //Absorber

  G4Box* solidA = new G4Box("Abso",worldX,worldX,absLength*0.5);
  logicA2 = new G4LogicalVolume( solidA,absMaterial,"Abs2");
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,absZ2),
			 "Abs2",logicA2,world,false,0);

  G4cout << "Absorber is " << G4BestUnit(absLength,"Length")
	 << " of " << absMaterial->GetName() << G4endl;

  //York

  G4Box* solidYV = new G4Box("VolY",worldX,worldX,york*0.5+absLength);
  logicYV = new G4LogicalVolume( solidYV,yorkMaterial,"VolY");
  G4VPhysicalVolume* physYV = new G4PVPlacement(0,G4ThreeVector(0.,0.,yorkZ),
						"VolY",logicYV,world,false,0);

  G4Box* solidY = new G4Box("York",worldX,worldX,york*0.5);
  logicY = new G4LogicalVolume( solidY,yorkMaterial,"York");
  pv = new G4PVPlacement(0,G4ThreeVector(),
			 "York",logicY,physYV,false,0);

  logicA3 = new G4LogicalVolume( solidA,absMaterial,"Abs3");
  logicA4 = new G4LogicalVolume( solidA,absMaterial,"Abs4");

  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,-(york+absLength)*0.5),
			 "Abs3",logicA3,physYV,false,0);
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,(york+absLength)*0.5),
			 "Abs4",logicA4,physYV,false,0);

  //Vertex volume
  G4Box* solidVV = new G4Box("VolV",worldX,worldX,vertexLength*2.+absLength+gap);
  logicVV = new G4LogicalVolume( solidVV,worldMaterial,"VolV");
  G4VPhysicalVolume* physVV = new G4PVPlacement(0,G4ThreeVector(0.,0.,vertexZ),
						"VolV",logicVV,world,false,0);

  //Absorber
  logicA1 = new G4LogicalVolume( solidA,absMaterial,"Abs1");
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,vertexLength*2.-absLength*0.5),
			 "Abs1",logicA1,physVV,false,0);

  //Vertex
  G4double vertWidth = ecalWidth/5.;
  G4int npads = (G4int)(vertWidth/padWidth);
  //G4cout << " vertWidth= " << vertWidth << " padWidth= " << padWidth 
  //	 << " npads= " << npads << G4endl;
  // insure beam to hit a middle of central pad
  npads = (npads/2)*2 + 1;
  x0 = -0.5*(padWidth + vertWidth);
  G4double x1 = 0.5*vertWidth + gap; 
  G4double z  = -(vertexLength+absLength);

  G4Box* solidVD = new G4Box("VertDet",x1,ecalWidth*0.5+gap,padLength*0.5);
  logicVD = new G4LogicalVolume( solidVD,vertMaterial,"VertDet");
  logicVD->SetSolid(solidVD);

  G4Box* solidV = new G4Box("Vert",padWidth*0.5,ecalWidth*0.5,padLength*0.5);
  logicV = new G4LogicalVolume( solidV,vertMaterial,"Vert");

  for (i=0; i<3; i++) {
    pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,z),"VertDet",logicVD,
                                    physVV,false,i);
    z += vertexLength;
  }
  x = x0;

  for (j=0; j<npads; j++) {

    pv = new G4PVPlacement(0,G4ThreeVector(x,0.,0.),logicV,"Vert",logicVD,
                                    false,k);
    x += padWidth;
  }

  G4cout << "Vertex is " << G4BestUnit(vertexLength,"Length")
         << " of 3 layers of Si of " << G4BestUnit(padLength,"Length")
         << " npads= " << npads
	 << G4endl;

  // Define region for the vertex detector
  vertexRegion->AddRootLogicalVolume(logicVV);
  vertexRegion->AddRootLogicalVolume(logicA3);

  // Define region for the muon detector
  muonRegion->AddRootLogicalVolume(logicYV);

  // color regions
  logicVV-> SetVisAttributes(G4VisAttributes::Invisible);
  logicV-> SetVisAttributes(G4VisAttributes::Invisible);
  logicECal-> SetVisAttributes(G4VisAttributes::Invisible);
  logicYV-> SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* regWcolor = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  logicWorld->SetVisAttributes(regWcolor);

  G4VisAttributes* regVcolor = new G4VisAttributes(G4Colour(0., 0.3, 0.7));
  logicVD->SetVisAttributes(regVcolor);

  G4VisAttributes* regCcolor = new G4VisAttributes(G4Colour(0., 0.7, 0.3));
  logicCal->SetVisAttributes(regCcolor);

  G4VisAttributes* regAcolor = new G4VisAttributes(G4Colour(1., 0.5, 0.5));
  logicA1->SetVisAttributes(regAcolor);
  logicA2->SetVisAttributes(regAcolor);
  logicA3->SetVisAttributes(regAcolor);
  logicA4->SetVisAttributes(regAcolor);

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
    G4NistManager::Instance()->FindOrBuildMaterial(mat, false);
  if (pttoMaterial) {
    calMaterial = pttoMaterial;
    if(logicCal) {
      logicCal->SetMaterial(calMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* pttoMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial(mat, false);
  if (pttoMaterial) {
    absMaterial = pttoMaterial;
    if(logicA1) {
      logicA1->SetMaterial(absMaterial);
      logicA2->SetMaterial(absMaterial);
      logicA3->SetMaterial(absMaterial);
      logicA4->SetMaterial(absMaterial);
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
    ecalLength = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetEcalWidth  (G4double val)   
{
  if(val > 0.0) {
    ecalWidth = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetVertexLength (G4double val) 
{
  if(val > 0.0) {
    vertexLength = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPadLength  (G4double val)   
{
  if(val > 0.0) {
    padLength = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPadWidth  (G4double val)    
{
  if(val > 0.0) {
    padWidth = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsLength(G4double val)     
{
  if(val > 0.0) {
    absLength = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
