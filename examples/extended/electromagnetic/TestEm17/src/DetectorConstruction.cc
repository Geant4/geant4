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
/// \file electromagnetic/TestEm17/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),fPBox(0), fLBox(0), fMaterial(0),
 fDetectorMessenger(0)
{
  fBoxSize = 1*m;
  DefineMaterials();
  SetMaterial("G4_Fe");  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  //
  // define Elements
  //
  G4double z,a;
  
  G4Element* C  = new G4Element("Carbon"   ,"C" , z= 6., a=  12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   ,"O" , z= 8., a=  16.00*g/mole);
  G4Element* Ca = new G4Element("Calcium"  ,"Ca", z=20., a=  40.08*g/mole);
  G4Element* H  = new G4Element("Hydrogen"   ,"H" , z= 1., a=  1.01*g/mole);
  G4Element* I  = new G4Element("Iodine"   ,"I" , z= 53., a=  126.90*g/mole);
  
  //
  // define materials
  //
  G4double density;
  G4int ncomponents, natoms;
  G4double fractionmass;
    
  new G4Material("galactic", z=1., a= 1.01*g/mole, universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);
    
  G4Material* Air = 
  new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=70.*perCent);
  Air->AddElement(O, fractionmass=30.*perCent);

  G4Material* CaCO3 = 
  new G4Material("CaCO3", density= 2.80*g/cm3, ncomponents=3);
  CaCO3->AddElement(Ca, natoms= 1);
  CaCO3->AddElement(C , natoms= 1);
  CaCO3->AddElement(O , natoms= 3);
   
  new G4Material("Carbon", z= 6., a= 12.01*g/mole, density= 2.265*g/cm3);   
  new G4Material("Iron"  , z=26., a= 55.85*g/mole, density= 7.870*g/cm3);
  new G4Material("Tin"   , z=50., a= 118.7*g/mole, density= 7.310*g/cm3);

  G4Material* HI = 
  new G4Material("HI", density= 2.8*g/cm3, ncomponents=2); //-35oC
  HI->AddElement(H, natoms= 1);
  HI->AddElement(I, natoms= 1);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  if(fPBox) { return fPBox; }

  G4Box*
  sBox = new G4Box("Container",                           //its name
                   fBoxSize/2,fBoxSize/2,fBoxSize/2);     //its dimensions

  fLBox = new G4LogicalVolume(sBox,                       //its shape
                             fMaterial,                   //its material
                             fMaterial->GetName());       //its name

  fPBox = new G4PVPlacement(0,                            //no rotation
                             G4ThreeVector(),             //at (0,0,0)
                           fLBox,                         //its logical volume
                           fMaterial->GetName(),          //its name
                           0,                             //its mother  volume
                           false,                         //no boolean operation
                           0);                            //copy number
                           
  PrintParameters();
  
  //always return the root volume
  //
  return fPBox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box is " << G4BestUnit(fBoxSize,"Length")
         << " of " << fMaterial->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if(pttoMaterial &&  fMaterial != pttoMaterial) {
    fMaterial = pttoMaterial;
    if(fLBox) { fLBox->SetMaterial(fMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  fBoxSize = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
