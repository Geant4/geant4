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
/// \file exoticphysics/monopole/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
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
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4UserLimits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4StateManager.hh"

#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4MonopoleFieldSetup.hh"
//#include "G4FieldManager.hh"
//#include "G4TransportationManager.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh" 
#include "G4SystemOfUnits.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fWorldMaterial(0),           
   fAbsorMaterial(0),
   fLogAbsor(0),
   fMonFieldSetup(),
   fDetectorMessenger(0)
{
  // default parameter values
  fAbsorSizeX = fAbsorSizeYZ = 10 * cm;
  fWorldSizeX = fWorldSizeYZ = 1.2 * fAbsorSizeX;
  fMaxStepSize = 5 * mm;

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
  //  delete fMonFieldSetup;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  /****************************    World   *****************************/
  G4Box * sWorld = new G4Box("world",                                
               fWorldSizeX / 2, fWorldSizeYZ / 2, fWorldSizeYZ / 2);

  G4LogicalVolume * lWorld = new G4LogicalVolume(sWorld,
                                                 fWorldMaterial,
                                                 "world");        

  G4VPhysicalVolume * pWorld = new G4PVPlacement(0,      //no rotation
                             G4ThreeVector(),            //at (0,0,0)
                           lWorld,                       //logical volume
                           "world",                      //name
                           0,                            //mother  volume
                           false,                        //no boolean operation
                           0);                           //copy number


  /**************************    Absorber    ***************************/
  G4Box * sAbsor = new G4Box("Absorber",                   
          fAbsorSizeX / 2, fAbsorSizeYZ / 2, fAbsorSizeYZ / 2);        

  fLogAbsor = new G4LogicalVolume(sAbsor,        
                                  fAbsorMaterial,
                                  "Absorber");
  
  new G4PVPlacement(0,                                //no rotation
                    G4ThreeVector(),                //at (0,0,0)
                    fLogAbsor,                        //logical volume
                    "Absorber",                        //name
                    lWorld,                               //mother  volume
                    false,                        //no boolean operation
                    0);                                //copy number
  fLogAbsor->SetUserLimits(new G4UserLimits(fMaxStepSize));

  PrintParameters();

  /************     always return the World volume     *****************/
  return pWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The Absorber is " << G4BestUnit(fAbsorSizeX, "Length")
         << " of " << fAbsorMaterial->GetName() << G4endl;
  G4cout << "\n---------------------------------------------------------\n";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeX(G4double value)
{
  if(value > 0.0) {
    fAbsorSizeX = value; 
    fWorldSizeX = 1.2 * fAbsorSizeX;
    if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
      G4RunManager::GetRunManager()->ReinitializeGeometry();
    }
  }
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeYZ(G4double value)
{
  if(value > 0.0) {
    fAbsorSizeYZ = value; 
    fWorldSizeYZ = 1.2 * fAbsorSizeYZ;
    if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
      G4RunManager::GetRunManager()->ReinitializeGeometry();
    }
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::ConstructSDandField()
{
  // Define magnetic field
  if ( ! fMonFieldSetup.Get() ) {
    G4MonopoleFieldSetup* fieldSetup
      = new G4MonopoleFieldSetup();
    G4AutoDelete::Register(fieldSetup); // Kernel will delete the F01FieldSetup
    fMonFieldSetup.Put(fieldSetup);
  }
  fMonFieldSetup.Get()->ConstructMagField(); // add field value
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaxStepSize(G4double step)
{
  fMaxStepSize = step;
  if(fLogAbsor) { fLogAbsor->SetUserLimits(new G4UserLimits(fMaxStepSize)); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
