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
/// \file medical/fanoCavity2/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : fWall(nullptr), fCavity(nullptr)
{
  // default parameter values
  fWallThickness   = 5*cm;  
  fCavityThickness = 2*mm;
  fCavityRadius    = 10*m;      
   
  DefineMaterials();
  SetWallMaterial("G4_WATER");
  SetCavityMaterial("g4Water_gas");
    
  // create commands for interactive definition of the detector  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4double z,a;
  
  G4Element* H  = new G4Element("Hydrogen" ,"H" , z= 1., a=   1.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   ,"O" , z= 8., a=  16.00*g/mole);

  G4Material* H2O = 
  new G4Material("Water", 1.0*g/cm3, 2);
  H2O->AddElement(H, 2);
  H2O->AddElement(O, 1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  
  G4Material* gas = 
  new G4Material("Water_gas", 1.0*mg/cm3, 2);
  gas->AddElement(H, 2);
  gas->AddElement(O, 1);
  gas->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  
  new G4Material("Graphite",     6, 12.01*g/mole, 2.265*g/cm3);
  new G4Material("Graphite_gas", 6, 12.01*g/mole, 2.265*mg/cm3);  
  
  new G4Material("Aluminium",     13, 26.98*g/mole, 2.700*g/cm3);
  new G4Material("Aluminium_gas", 13, 26.98*g/mole, 2.700*mg/cm3);
  
  // alternatively, use G4 data base
  //
  G4NistManager* nist = G4NistManager::Instance();

  nist->FindOrBuildMaterial("G4_WATER");
  nist->BuildMaterialWithNewDensity("g4Water_gas","G4_WATER", 1.0*mg/cm3);  
 
  // printout 
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if(fWall) { return fWall; }

  // Chamber
  //
  fTotalThickness = fCavityThickness + 2*fWallThickness;
  fWallRadius     = fCavityRadius;
  
  G4Tubs* 
  sChamber = new G4Tubs("Chamber",                                     //name
                         0.,fWallRadius,0.5*fTotalThickness,0.,twopi); //size

  G4LogicalVolume*
  lChamber = new G4LogicalVolume(sChamber,              //solid
                                    fWallMaterial,     //material
                                    "Chamber");         //name
                                   
  fWall = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(),            //at (0,0,0)
                            lChamber,                   //logical volume
                            "Wall",                     //name
                            0,                          //mother  volume
                            false,                      //no boolean operation
                            0);                         //copy number

  // Cavity
  //                             
  G4Tubs*
  sCavity = new G4Tubs("Cavity",        
                       0.,fCavityRadius,0.5*fCavityThickness,0.,twopi);
                 
  G4LogicalVolume* 
  lCavity = new G4LogicalVolume(sCavity,                //shape
                                fCavityMaterial,        //material
                                "Cavity");              //name
                                
  fCavity = new G4PVPlacement(0,                        //no rotation
                             G4ThreeVector(),           //at (0,0,0)
                             lCavity,                   //logical volume
                            "Cavity",                   //name
                             lChamber,                  //mother  volume
                             false,                     //no boolean operation
                             1);                        //copy number
                                
  PrintParameters();
    
  //
  //always return the root volume
  //  
  return fWall;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Wall is " << G4BestUnit(fWallThickness,"Length")
         << " of " << fWallMaterial->GetName() << " ( " 
         << G4BestUnit(fWallMaterial->GetDensity(),"Volumic Mass") << " )\n";
  G4cout << "     The Cavity is " << G4BestUnit(fCavityThickness,"Length")
         << " of " << fCavityMaterial->GetName() << " ( " 
         << G4BestUnit(fCavityMaterial->GetDensity(),"Volumic Mass") << " )";
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWallThickness(G4double value)
{
  fWallThickness = value;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWallMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);         
  if (pttoMaterial) fWallMaterial = pttoMaterial;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCavityThickness(G4double value)
{
  fCavityThickness = value;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCavityRadius(G4double value)
{
  fCavityRadius  = value;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCavityMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);         
  if (pttoMaterial) fCavityMaterial = pttoMaterial; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
