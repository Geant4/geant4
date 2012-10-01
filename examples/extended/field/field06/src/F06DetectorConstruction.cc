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
/// \file field/field06/src/F06DetectorConstruction.cc
/// \brief Implementation of the F06DetectorConstruction class
//
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F06DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"

#include "F06Field.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F06DetectorConstruction::F06DetectorConstruction()
 : Vacuum(0), field(0)
{
  // materials
  DefineMaterials();

  // ensure the global field is initialized
  field = new F06Field();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F06DetectorConstruction::~F06DetectorConstruction()
{
  delete field;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F06DetectorConstruction::DefineMaterials()
{
  G4NistManager* nistMan = G4NistManager::Instance();

  Vacuum = nistMan->FindOrBuildMaterial("G4_Galactic");

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* F06DetectorConstruction::Construct()
{
  //     
  // World
  //

  G4double expHall_x = 1.0*m;
  G4double expHall_y = 1.0*m;
  G4double expHall_z = 1.0*m;

  solidWorld = new G4Box("World",                  //its name
                   expHall_x,expHall_y,expHall_z); //its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,     //its solid
                                   Vacuum,         //its material
                                   "World");       //its name
                                   
  physiWorld = new G4PVPlacement(0,                //no rotation
                                 G4ThreeVector(),  //at (0,0,0)
                                 logicWorld,       //its logical volume
                                 "World",          //its name
                                 0,                //its mother  volume
                                 false,            //no boolean operation
                                 0);               //copy number
  
  G4double maxStep = 1.0*mm;
  G4double maxTime = 41.*s;

  G4UserLimits* stepLimit = new G4UserLimits(maxStep,DBL_MAX,maxTime);

  logicWorld->SetUserLimits(stepLimit);
 
  //                                        
  // Visualization attributes
  //
  // logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

  //
  //always return the physical World
  //
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
