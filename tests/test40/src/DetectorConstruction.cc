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
// $Id$

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4NistManager.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4UniformMagField* DetectorConstruction::magField = 0;

DetectorConstruction::DetectorConstruction()
:nLtot(20),nRtot(20),dLradl(1.),dRradl(0.25),
 myMaterial(0),
 EcalLength(0.),EcalRadius(0.)    ,
 solidEcal(0) ,logicEcal(0) ,physiEcal(0),
 solidSlice(0),logicSlice(0),physiSlice(0),
 solidRing(0) ,logicRing(0) ,physiRing(0)
{
  DefineMaterials();
  detectorMessenger = new DetectorMessenger(this);
  edeptrue = 0.9425;
  rmstrue = 0.0355;
  limittrue = DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete detectorMessenger;
  delete magField;
  magField = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  //
  // define few Elements
  //
  G4double a, z;
    
  G4Element* H  = new G4Element("Hydrogen",  "H", z= 1., a=   1.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",  "N", z= 7., a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,  "O", z= 8., a=  16.00*g/mole);
  G4Element* Ge = new G4Element("Germanium", "Ge",z=32., a=  72.59*g/mole);
  G4Element* W  = new G4Element("Tungsten",  "W", z=74., a= 183.84*g/mole);
  G4Element* Pb = new G4Element("Lead",      "Pb",z=82., a= 207.19*g/mole);
  G4Element* Bi = new G4Element("Bismuth",   "Bi",z=83., a= 208.98*g/mole);

  //
  // define materials
  //
  G4double density;
  G4double fractionmass;  G4int ncomponents, natoms;

  G4Material* Air = 
  new G4Material("Air", density= 1.29*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  G4Material* H2O = 
  new G4Material("Water", density= 1.00*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);

  new G4Material("Aluminium",   z=13., a= 26.98*g/mole, density= 2.7*g/cm3);

  G4Material* Fe = 
  new G4Material("Iron",        z=26., a= 55.85*g/mole, density= 7.87*g/cm3);
  G4Material* Ni = 
  new G4Material("Nickel",      z=28., a= 58.69*g/mole, density= 8.96*g/cm3);
  G4Material* Cu = 
  new G4Material("Copper",      z=29., a= 63.54*g/mole, density= 8.96*g/cm3);

  G4Material* BGO = 
  new G4Material("BGO", density= 7.10*g/cm3, ncomponents=3);
  BGO->AddElement(O , natoms=12);
  BGO->AddElement(Ge, natoms= 3);
  BGO->AddElement(Bi, natoms= 4);

  G4Material* PbWO = 
  new G4Material("PbWO4", density= 8.28*g/cm3, ncomponents=3);
  PbWO->AddElement(O , natoms=4);
  PbWO->AddElement(Pb, natoms=1);
  PbWO->AddElement(W , natoms=1);

  G4Material* w = 
  new G4Material("Tungsten", density= 19.30*g/cm3, ncomponents=1);
  w->AddElement(W, fractionmass=1.0);

  G4Material* ma1 = new G4Material("FCal2Slugs",density = 18.6*g/cm3, 3);
  ma1->AddMaterial(w, fractionmass=0.97);
  ma1->AddMaterial(Fe,fractionmass=0.01);
  ma1->AddMaterial(Ni,fractionmass=0.02);

  G4Material* ma2 = new G4Material("FCal2Abs",density = 10.*g/cm3, 2);
  ma2->AddMaterial(Cu, fractionmass=0.2);
  ma2->AddMaterial(ma1,fractionmass=0.8);

  G4Material* pb = 
  new G4Material("Lead", density= 11.35*g/cm3, ncomponents=1);
  pb->AddElement(Pb, fractionmass=1.0);

  myMaterial    = PbWO;
  worldMaterial = Air;

  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  G4double Radl = myMaterial->GetRadlen();

  G4double dL = dLradl*Radl, dR = dRradl*Radl;  
  EcalLength = nLtot*dL; EcalRadius = nRtot*dR;

  //
  // World
  //
  G4Box* solidW = new G4Box("World",EcalRadius,EcalRadius,EcalLength);
  G4LogicalVolume* logicW = new G4LogicalVolume( solidW,worldMaterial,
                                                "World");
  G4VPhysicalVolume* world = new G4PVPlacement(0,G4ThreeVector(),
                                       "World",logicW,0,false,0);
  
  //Vertex volume

  G4Box* solidVV = 
    new G4Box("VolV",EcalRadius*0.8,EcalRadius*0.8,EcalLength*0.01);
  G4LogicalVolume* logicVV = new G4LogicalVolume( solidVV,worldMaterial,"VolV");
  G4VPhysicalVolume* physVV = 
    new G4PVPlacement(0,G4ThreeVector(0.,0.,-EcalLength*0.6),
		      "VolV",logicVV,world,false,0);


  // Vertex detector

  G4int    nl = 1000;
  G4double xl = (G4double)nl;
  G4double widthX = 2.0*EcalRadius*0.8/xl;
  G4double x0 = -0.5*widthX*(xl - 1.0);

  G4Box* solidLayer = 
    new G4Box("Layer",widthX*0.49,EcalRadius*0.8,EcalLength*0.0099);
  G4LogicalVolume* logicLayer = 
    new G4LogicalVolume(solidLayer,worldMaterial,"Layer");       
  for (G4int j=0; j<nl; j++) {
    new G4PVPlacement(0,G4ThreeVector(x0,0.,0.),"Layer",
		      logicLayer,physVV,false,j);
    x0 += widthX;
  }
  
  //   
  // Ecal
  //
  solidEcal = new G4Tubs("Ecal",0.,EcalRadius,0.5*EcalLength,0.,360*deg);
  logicEcal = new G4LogicalVolume( solidEcal,myMaterial,"Ecal",0,0,0);
  physiEcal = 
    new G4PVPlacement(0,G4ThreeVector(),"Ecal",logicEcal,world,false,0);
				
  // Ring
  //
  for (G4int i=0; i<nRtot; i++)
     {
      solidRing = new G4Tubs("Ring",i*dR,(i+1)*dR,0.5*EcalLength,0.,360*deg);
      logicRing = new G4LogicalVolume(solidRing,myMaterial,"Ring",0,0,0);
      physiRing = new G4PVPlacement(0,G4ThreeVector(),"Ring",logicRing,
                                    physiEcal,false,i);
				                   				
      // Slice
      solidSlice = new G4Tubs("Slice",i*dR,(i+1)*dR,0.5*dL,0.,360*deg);
      logicSlice = new G4LogicalVolume(solidSlice,myMaterial,"Slice",0,0,0);
      if (nLtot >1)
        physiSlice = new G4PVReplica("Slice",logicSlice,physiRing,
                                    kZAxis,nLtot,dL);
      else
        physiSlice = new G4PVPlacement(0,G4ThreeVector(),"Slice",logicSlice,
                                    physiRing,false,0);
     }				                   

  G4cout << "Absorber is " << G4BestUnit(EcalLength,"Length") 
       << " of " << myMaterial->GetName() << G4endl; 

  //Muon detector volume

  G4Box* solidM = new G4Box("VolM",EcalRadius*0.8,EcalRadius*0.8,EcalLength*0.1);
  G4LogicalVolume* logicM = new G4LogicalVolume( solidM,myMaterial,"VolM");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,EcalLength*0.8),
                                       "VolM",logicM,world,false,0);

  //
  //always return the physical World
  //
  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if(pttoMaterial && myMaterial != pttoMaterial) {
    myMaterial = pttoMaterial;
    if(logicEcal)  { logicEcal->SetMaterial(myMaterial); }
    if(logicRing)  { logicRing->SetMaterial(myMaterial); }
    if(logicSlice) { logicSlice->SetMaterial(myMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& materialChoice)
{
    // search the material by its name
    G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
    if (pttoMaterial) {
        worldMaterial = pttoMaterial;
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetLBining(G4ThreeVector Value)
{
  nLtot = (G4int)Value(0);
  dLradl = Value(1);
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetRBining(G4ThreeVector Value)
{
  nRtot = (G4int)Value(0);
  dRradl = Value(1);
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetEdepAndRMS(G4ThreeVector Value)
{
  edeptrue = Value(0);
  rmstrue  = Value(1);
  limittrue= Value(2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    if ( fFieldMessenger.Get() == 0 ) {
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue = G4ThreeVector();
        G4GlobalMagFieldMessenger* msg =
        new G4GlobalMagFieldMessenger(fieldValue);
        //msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );
        
    }
//
//  if(!magField && 0.0 != fieldValue) {
//    G4FieldManager* fieldMgr
//      = G4TransportationManager::GetTransportationManager()->GetFieldManager();
//    magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
//    fieldMgr->SetDetectorField(magField);
//    fieldMgr->CreateChordFinder(magField);
//  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
