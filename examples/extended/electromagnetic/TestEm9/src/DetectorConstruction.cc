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
// $Id: DetectorConstruction.cc,v 1.6 2004/12/02 19:06:05 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
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
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:logicC(0),
 logicA1(0),
 logicA2(0)
{
  detectorMessenger = new DetectorMessenger(this);
  DefineMaterials();
  ecalLength   = 36.*cm;
  ecalWidth    = 6.*cm;
  vertexLength = 3.*cm;
  padLength    = 0.1*mm;
  padWidth     = 0.02*mm;
  absLength    = 2.*mm;
  vertexRegion = 0;
  muonRegion   = 0;
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
  G4String name, symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4double fractionmass;

  //
  // define few Elements
  //

    a = 1.01*g/mole;
    G4Element* H = new G4Element(name="Hydrogen", symbol="H", z=1., a);

    a = 14.01*g/mole;
    G4Element* N = new G4Element(name="Nitrogen", symbol="N", z=7., a);

    a = 16.00*g/mole;
    G4Element* O = new G4Element(name="Oxygen"  , symbol="O", z=8., a);

    a = 72.59*g/mole;
    G4Element* Ge = new G4Element(name="Germanium", symbol="Ge",z=32., a);

    a = 183.84*g/mole;
    G4Element* W = new G4Element(name="Tungsten"  , symbol="W", z=74., a);

    a = 207.19*g/mole;
    G4Element* Pb = new G4Element(name="Lead"     , symbol="Pb",z=82., a);

    a = 208.98*g/mole;
    G4Element* Bi = new G4Element(name="Bismuth"  , symbol="Bi",z=83., a);

    G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);

    G4Element*   I  = new G4Element ("Iodide"  , "I", 53. , 126.9044*g/mole);

  //
  // define materials
  //

    G4Material* ma = 0;

  //Air
    density = 1.29*mg/cm3;
    G4Material* Air = new G4Material(name="Air", density, ncomponents=2);
    Air->AddElement(N, fractionmass=0.7);
    Air->AddElement(O, fractionmass=0.3);

  //H2O
    density = 1.00*g/cm3;
    G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
    H2O->AddElement(H, natoms=2);
    H2O->AddElement(O, natoms=1);
    G4double exc = H2O->GetIonisation()->FindMeanExcitationEnergy("H_2O");
    H2O->GetIonisation()->SetMeanExcitationEnergy(exc);

  //liquid argon
    a = 39.95*g/mole;
    density = 1.390*g/cm3;
    ma = new G4Material(name="lAr", z=18., a, density);

  //Al
    a = 26.98*g/mole;
    density = 2.7*g/cm3;
    absMaterial = new G4Material(name="Al", z=13., a, density);

    //Si
    density = 2.330*g/cm3;
    a = 28.09*g/mole;
    vertMaterial = new G4Material(name="Si", z=14., a, density);


    //Fe
    a = 55.85*g/mole;
    density = 7.87*g/cm3;
    yorkMaterial = new G4Material(name="Fe", z=26., a, density);

    // CsI
    ma = new G4Material ("CsI" , 4.51*g/cm3, 2);
    ma->SetChemicalFormula("CsI");
    ma->AddElement(Cs,1);
    ma->AddElement(I,1);
    ma->GetIonisation()->SetMeanExcitationEnergy(415.689*eV);
    calMaterial = ma;


  //BGO
    density = 7.10*g/cm3;
    G4Material* BGO = new G4Material(name="BGO", density, ncomponents=3);
    BGO->AddElement(O , natoms=12);
    BGO->AddElement(Ge, natoms= 3);
    BGO->AddElement(Bi, natoms= 4);

  //PbWO4
    density = 8.28*g/cm3;
    G4Material* PbWO = new G4Material(name="PbWO4", density, ncomponents=3);
    PbWO->AddElement(O , natoms=4);
    PbWO->AddElement(Pb, natoms=1);
    PbWO->AddElement(W , natoms=1);

  //Tungsten
    density = 19.30*g/cm3;
    G4Material* w = new G4Material(name="Tungsten", density, ncomponents=1);
    w->AddElement(W, fractionmass=1.0);

  //Pb
    density = 11.35*g/cm3;
    G4Material* pb = new G4Material(name="Lead", density, ncomponents=1);
    pb->AddElement(Pb, fractionmass=1.0);

    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //choose material
  worldMaterial = Air;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  if(vertexRegion) delete vertexRegion;
  if(muonRegion) delete muonRegion;
  vertexRegion = new G4Region("VertexDetector");
  muonRegion   = new G4Region("MuonDetector");

  if(vertexLength < padLength*5.0) vertexLength = padLength*5.0;
  G4double gap    = 0.01*mm;
  G4double biggap = 2.*cm;
  G4double york   = 10.*cm;

           worldZ = 2.*vertexLength + 3.*absLength + 0.5*(ecalLength + york) + biggap*2.;
  G4double worldX = ecalWidth*3.0;
  G4double vertexZ= -worldZ + vertexLength*2.0 + absLength     + biggap;
  G4double absZ2  = -worldZ + vertexLength*4.0 + absLength*3.5 + biggap;
  G4double ecalZ  = -worldZ + vertexLength*4.0 + absLength*4.0 + ecalLength*0.5 + 2.*biggap;
  G4double yorkZ  = -worldZ + vertexLength*4.0 + absLength*5.0 + ecalLength
                            + york*0.5 + 3.*biggap;

  //
  // World
  //
  G4Box* solidW = new G4Box("World",worldX,worldX,worldZ);
  G4LogicalVolume* logicW = new G4LogicalVolume( solidW,worldMaterial,
                                                "World");
  G4VPhysicalVolume* world = new G4PVPlacement(0,G4ThreeVector(),
                                       "World",logicW,0,false,0);

  //
  // Ecal
  //
  G4Box* solidE = new G4Box("VolE",worldX,worldX,ecalLength*0.5 + gap);
  G4LogicalVolume* logicE = new G4LogicalVolume( solidE,worldMaterial,
                                                "VolE");
  G4VPhysicalVolume* physE = new G4PVPlacement(0,G4ThreeVector(0.,0.,ecalZ),
                                       "VolE",logicE,world,false,0);

  G4Box* solidC = new G4Box("Ecal",ecalWidth*0.5,ecalWidth*0.5,ecalLength*0.5);
  logicC = new G4LogicalVolume( solidC,calMaterial,"Ecal");

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

      pv = new G4PVPlacement(0,G4ThreeVector(x,y,0.),"Ecal",logicC,
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
  G4LogicalVolume* logicYV = new G4LogicalVolume( solidYV,yorkMaterial,"VolY");
  G4VPhysicalVolume* physYV = new G4PVPlacement(0,G4ThreeVector(0.,0.,yorkZ),
                                       "VolY",logicYV,world,false,0);

  G4Box* solidY = new G4Box("York",worldX,worldX,york*0.5);
  G4LogicalVolume* logicY = new G4LogicalVolume( solidY,yorkMaterial,"York");
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
  G4LogicalVolume* logicVV = new G4LogicalVolume( solidVV,worldMaterial,"VolV");
  G4VPhysicalVolume* physVV = new G4PVPlacement(0,G4ThreeVector(0.,0.,vertexZ),
                                       "VolV",logicVV,world,false,0);

  //Absorber

  logicA1 = new G4LogicalVolume( solidA,absMaterial,"Abs1");
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,vertexLength*2.-absLength*0.5),
                                       "Abs1",logicA1,physVV,false,0);

  //Vertex

  G4double vertWidth = ecalWidth/5.;
  G4int npads = (G4int)(vertWidth/padWidth);
  npads = (npads/2)*2 + 1;
  x0 = -0.5*padWidth*((G4double)(npads-1));
  G4double x1 = std::fabs(x0) + 0.5*padWidth + gap; 
  G4double z  = -(vertexLength+absLength);

  G4Box* solidVD = new G4Box("VertDet",x1,ecalWidth*0.5+gap,padLength*0.5);
  G4LogicalVolume* logicVD = new G4LogicalVolume( solidVD,vertMaterial,"VertDet");

  G4Box* solidV = new G4Box("Vert",padWidth*0.5,ecalWidth*0.5,padLength*0.5);
  G4LogicalVolume* logicV = new G4LogicalVolume( solidV,vertMaterial,"Vert");

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
  logicE-> SetVisAttributes(G4VisAttributes::Invisible);
  logicYV-> SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* regWcolor = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  logicW->SetVisAttributes(regWcolor);

  G4VisAttributes* regVcolor = new G4VisAttributes(G4Colour(0., 0.3, 0.7));
  logicVD->SetVisAttributes(regVcolor);

  G4VisAttributes* regCcolor = new G4VisAttributes(G4Colour(0., 0.7, 0.3));
  logicC->SetVisAttributes(regCcolor);

  G4VisAttributes* regAcolor = new G4VisAttributes(G4Colour(1., 0.5, 0.5));
  logicA1->SetVisAttributes(regAcolor);
  logicA2->SetVisAttributes(regAcolor);
  logicA3->SetVisAttributes(regAcolor);
  logicA4->SetVisAttributes(regAcolor);

  G4VisAttributes* regMcolor = new G4VisAttributes(G4Colour(1., 1., 0.));
  logicY->SetVisAttributes(regMcolor);

  // always return world
  //
  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetEcalMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial)
     {
        calMaterial = pttoMaterial;
        if(logicC) logicC->SetMaterial(calMaterial);
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial)
     {
        absMaterial = pttoMaterial;
        if(logicA1) logicA1->SetMaterial(absMaterial);
        if(logicA2) logicA2->SetMaterial(absMaterial);
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
