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
// $Id: Tst26DetectorConstruction.cc,v 1.3 2003-02-06 11:53:27 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst26DetectorConstruction.hh"
#include "Tst26DetectorMessenger.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4TransportationManager.hh"

#include "G4RunManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26DetectorConstruction::Tst26DetectorConstruction()
:ecalLength(20.*cm),
 ecalWidth(4.*cm),
 vertexLength(3.*cm),
 padLength(0.1*mm),
 padWidth(0.02*mm),
 absLength(2.*mm),
 logicC(0),
 logicA1(0),
 logicA2(0)
{
  detectorMessenger = new Tst26DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26DetectorConstruction::~Tst26DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Tst26DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26DetectorConstruction::DefineMaterials()
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
    //  G4double exc = H2O->GetIonisation()->FindMeanExcitationEnergy("H_2O");
    //  H2O->GetIonisation()->SetMeanExcitationEnergy(exc);
    
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
  calMaterial = PbWO;
  worldMaterial = Air;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* Tst26DetectorConstruction::ConstructVolumes()
{
  if(vertexLength < padLength*5.0) vertexLength = padLength*5.0;
  G4double gap    = 0.01*mm;
  G4double york   = 10*cm;

  G4double worldZ = vertexLength + absLength*10.5 + ecalLength*1.5 + york;
  G4double worldX = ecalWidth*3.0;
  G4double vertexZ= vertexLength*0.5 + absLength*2.5;
  G4double absZ2  = vertexLength     + absLength*6.0;
  G4double ecalZ  = vertexLength     + absLength*7.5 + ecalLength*0.5;
  G4double yorkZ  = vertexLength + absLength*9.5 + ecalLength + york*0.5;
  
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

  G4Box* solidVV = new G4Box("VolV",worldX,worldX,vertexLength*0.5+absLength*2.5);
  G4LogicalVolume* logicVV = new G4LogicalVolume( solidVV,worldMaterial,"VolV");
  G4VPhysicalVolume* physVV = new G4PVPlacement(0,G4ThreeVector(0.,0.,vertexZ),
                                       "VolV",logicVV,world,false,0);


  //Absorber

  logicA1 = new G4LogicalVolume( solidA,absMaterial,"Abs1");
  pv = new G4PVPlacement(0,G4ThreeVector(0.,0.,vertexLength*0.5+absLength*1.5),
                                       "Abs1",logicA1,physVV,false,0);

  //Vertex

  G4Box* solidV = new G4Box("Vert",padWidth*0.5,ecalWidth*0.5,padLength*0.5);
  G4LogicalVolume* logicV = new G4LogicalVolume( solidV,vertMaterial,"Vert");
  G4int npads = (G4int)(ecalWidth/padWidth);
  npads = (npads/2)*2 + 1;
  k = 0;
  x0 = -0.5*padWidth*((G4double)(npads-1));
  G4double z  = (vertexLength-padLength)*0.5;

  for (i=0; i<3; i++) {
    x = x0;
    y = z*(i - 1);
    for (j=0; j<npads; j++) {

      pv = new G4PVPlacement(0,G4ThreeVector(x,0.,y),"Vert",logicV,
                                    physVV,false,k);
      k++;
      x += padWidth;
    }
  }                  

  G4cout << "Vertex is " << G4BestUnit(vertexLength,"Length") 
         << " of 3 layers of Si of " << G4BestUnit(padLength,"Length") 
         << " npads= " << npads
    //         << " x0= " << x0
    //     << " z= " << z
	 << G4endl; 

  // Define region for the vertex detector

  G4RegionStore* rs = G4RegionStore::GetInstance();
  G4Region* regionV  = rs->GetRegion("VertexDetector");
  if( !regionV ) {
    regionV = new G4Region("VertexDetector");
    regionV->AddRootLogicalVolume(logicVV);
    regionV->AddRootLogicalVolume(logicA3);
  }

  // Define region for the muon detector

  G4Region* regionM  = rs->GetRegion("MuonDetector");
  if( !regionM ) {
    regionM = new G4Region("MuonDetector");
    regionM->AddRootLogicalVolume(logicYV);
  }

  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26DetectorConstruction::SetEcalMaterial(const G4String& materialChoice)
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

void Tst26DetectorConstruction::SetAbsMaterial(const G4String& materialChoice)
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
  
void Tst26DetectorConstruction::UpdateGeometry()
{
  G4VPhysicalVolume* v = ConstructVolumes();
  G4RunManager* rm = G4RunManager::GetRunManager();
  rm->GeometryHasBeenModified();
  rm->DefineWorldVolume(v);
  rm->ResetNavigator();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
