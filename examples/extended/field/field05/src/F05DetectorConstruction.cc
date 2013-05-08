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
// $Id$
//
/// \file field/field05/src/F05DetectorConstruction.cc
/// \brief Implementation of the F05DetectorConstruction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F05DetectorConstruction.hh"

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

#include "F05Field.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05DetectorConstruction::F05DetectorConstruction()
 : fVacuum(0), fWorldSizeXY(0), fWorldSizeZ(0), 
   fSolidWorld(0), fLogicWorld(0), fPhysiWorld(0), fField(0)
{
  // materials
  DefineMaterials();

  // ensure the global field is initialized
  fField = new F05Field();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05DetectorConstruction::~F05DetectorConstruction()
{
  delete fField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05DetectorConstruction::DefineMaterials()
{
  G4NistManager* nistMan = G4NistManager::Instance();

  fVacuum = nistMan->FindOrBuildMaterial("G4_Galactic");

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* F05DetectorConstruction::Construct()
{
  //
  // World
  //

  fWorldSizeXY = 20.0*m;
  fWorldSizeZ  =  1.0*mm;

  fSolidWorld = new G4Box("World",                               //its name
                   fWorldSizeXY/2,fWorldSizeXY/2,fWorldSizeZ/2); //its size
 
  fLogicWorld = new G4LogicalVolume(fSolidWorld,        //its solid
                                    fVacuum,            //its material
                                    "World");           //its name
 
  fPhysiWorld = new G4PVPlacement(0,                    //no rotation
                                  G4ThreeVector(),      //at (0,0,0)
                                  fLogicWorld,          //its logical volume
                                  "World",              //its name
                                  0,                    //its mother  volume
                                  false,                //no boolean operation
                                  0);                   //copy number
  
  G4UserLimits* stepLimit;
  stepLimit = new G4UserLimits(5*mm);

  fLogicWorld->SetUserLimits(stepLimit);
 
  //
  // Visualization attributes
  //
  // fLogicWorld->SetVisAttributes (G4VisAttributes::Invisible);

  //
  //always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
