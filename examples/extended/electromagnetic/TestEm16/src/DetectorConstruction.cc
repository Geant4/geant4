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
/// \file electromagnetic/TestEm16/src/DetectorConstruction.cc
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
#include "G4UserLimits.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fLBox(0), fBox(0), fBoxSize(500*m), fMaterial(0),
 fUserLimits(0), fDetectorMessenger(0)
{
  DefineMaterials();
  SetMaterial("Iron");

  // create UserLimits
  fUserLimits = new G4UserLimits();

  // create commands for interactive definition of the detector
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
 G4double a, z, density;

 new G4Material("Beryllium", z= 4., a= 9.012182*g/mole, density= 1.848*g/cm3);
 new G4Material("Carbon",    z= 6., a= 12.011*g/mole,   density= 2.265*g/cm3);
 new G4Material("Iron",      z=26., a= 55.85*g/mole,    density= 7.870*g/cm3);

 // define a vacuum with a restgas pressure  typical for accelerators
 G4double const Torr  = atmosphere/760.;         // 1 Torr
 G4double pressure = 10e-9*Torr, temperature = 296.150*kelvin;    // 23  Celsius
 new G4Material("Vacuum", z=7., a=14.01*g/mole, density= 1.516784e-11*kg/m3,
                 kStateGas, temperature, pressure);

 G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4Box*
  sBox = new G4Box("Container",                                //its name
                   fBoxSize/2,fBoxSize/2,fBoxSize/2);        //its dimensions

  fLBox = new G4LogicalVolume(sBox,                        //its shape
                             fMaterial,                        //its material
                             fMaterial->GetName());        //its name
                        
  fLBox->SetUserLimits(fUserLimits);                        

  fBox = new G4PVPlacement(0,                                //no rotation
                             G4ThreeVector(),                //at (0,0,0)
                           fLBox,                        //its logical volume
                           fMaterial->GetName(),        //its name
                           0,                            //its mother  volume
                           false,                        //no boolean operation
                           0);                                //copy number

  PrintParameters();

  //
  //always return the root volume
  //
  return fBox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box is " << G4BestUnit(fBoxSize,"Length")
         << " of " << fMaterial->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) {
     fMaterial = pttoMaterial;
     if ( fLBox ) fLBox->SetMaterial(fMaterial);
 }
 G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  fBoxSize = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaxStepSize(G4double val)
{
  // set the maximum allowed step size
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetMaxStepSize: maxStep "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  fUserLimits->SetMaxAllowedStep(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaxStepLength(G4double val)
{
  // set the maximum length of tracking step
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetMaxStepLength: maxStep "
             << val  << " out of range. Command refused" << G4endl;
      return;
    }
  G4TransportationManager* tmanager =
    G4TransportationManager::GetTransportationManager();
  tmanager->GetPropagatorInField()
          ->SetLargestAcceptableStep(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
