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
/// \file exoticphysics/dmparticle/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 68036 2013-03-13 14:13:45Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "TargetSD.hh"
#include "G4SDManager.hh"

#include "G4GeometryManager.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4FieldManager.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh" 
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorMaterial(nullptr),
   fLogAbsor(nullptr),
   fWorld(nullptr)
{
  // default parameter values
  fAbsorSizeZ = fAbsorSizeXY = 10 * cm;
  fWorldSizeZ = fWorldSizeXY = 1.2 * fAbsorSizeZ;

  SetMaterial("G4_Al");
  fWorldMaterial = 
    G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  // create commands for interactive definition of the detector
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  if(fWorld) { return fWorld; }

  /****************************    World   *****************************/
  G4Box * sWorld = new G4Box("world",                                
               fWorldSizeXY / 2, fWorldSizeXY / 2, fWorldSizeZ / 2);

  G4LogicalVolume * lWorld = new G4LogicalVolume(sWorld,
                                                 fWorldMaterial,
                                                 "world");        

  fWorld = new G4PVPlacement(0,               //no rotation
                             G4ThreeVector(), //at (0,0,0)
                             lWorld,          //logical volume
                             "world",         //name
                             0,               //mother  volume
                             false,           //no boolean operation
                             0);              //copy number


  /**************************    Absorber    ***************************/

  G4Box * sAbsor = new G4Box("Absorber",                   
          fAbsorSizeXY / 2, fAbsorSizeXY / 2, fAbsorSizeZ / 2);        

  fLogAbsor = new G4LogicalVolume(sAbsor,        
                                  fAbsorMaterial,
                                  "Absorber");
  
  new G4PVPlacement(0,               //no rotation
                    G4ThreeVector(), //at (0,0,0)
                    fLogAbsor,       //logical volume
                    "Absorber",      //name
                    lWorld,          //mother  volume
                    false,           //no boolean operation
                    0);              //copy number

  PrintParameters();

  /************     always return the World volume     *****************/
  return fWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Absorber is " << G4BestUnit(fAbsorSizeZ, "Length")
         << " of " << fAbsorMaterial->GetName() << G4endl;
  G4cout << "\n---------------------------------------------------------\n";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeZ(G4double value)
{
  if(value > 0.0) {
    fAbsorSizeZ = value; 
    fWorldSizeZ = 1.2 * fAbsorSizeZ;
  }
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeXY(G4double value)
{
  if(value > 0.0) 
  {
    fAbsorSizeXY = value; 
    fWorldSizeXY = 1.2 * fAbsorSizeXY;
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& namemat)
{
  // search the material by its name   

  G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial(namemat);
  if(!mat) {
    G4cout << "!!! DetectorConstruction::SetMaterial: WARNING Material <"
           << namemat << "> does not exist in DB" << G4endl;
    return;
  }
  // new material is found out
  if (mat != fAbsorMaterial) {
    fAbsorMaterial = mat;
    if(fLogAbsor) { fLogAbsor->SetMaterial(mat); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  auto sd = new TargetSD("Target");
  G4SDManager::GetSDMpointer()->AddNewDetector(sd);
  SetSensitiveDetector(fLogAbsor, sd); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
