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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction(),  
                                               fTargetMaterial( nullptr ),
  fLogicExperimentalHall( nullptr ), fPhysExperimentalHall( nullptr ),
  fLogicTargetLayer( nullptr ), fPhysTargetLayer( nullptr ),
  fFieldMgr( nullptr ), fUniformMagField( nullptr ),
  fDetectorMessenger( nullptr ),
  fTargetInnerRadius( 9.0*mm ), fTargetOuterRadius( 11.0*mm )  //***LOOKHERE*** Default values
{
  fFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  //***LOOKHERE*** Default material
  fTargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Be" );
  fDetectorMessenger = new DetectorMessenger( this );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
  delete fUniformMagField;
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {
  return ConstructLayer();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructLayer() {
  // Clean old geometry, if any
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  // The geometry consists of a cylinder with axis along the z-direction, and
  // positioned at the center, (0.0, 0.0, 0.0). Its inner and outer radius, and
  // its material can be set via UI commands.
  // The world volume (experimental hall) is a box slightly bigger than the cylinder
  // and it is filled of "G4_Galactic" material.
  const G4double halfLength = 1.0*m;  //***LOOKHERE*** Half-length of the cylinder
  const G4double expHall_x = 1.01*halfLength;  // half dimension along x
  const G4double expHall_y = 1.01*halfLength;  // half dimension along y
  const G4double expHall_z = 1.01*halfLength;  // half dimension along z
  G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Galactic" );
  G4Box* experimentalHallBox = new G4Box( "expHallBox", expHall_x, expHall_y, expHall_z );
  fLogicExperimentalHall = new G4LogicalVolume( experimentalHallBox,  // solid
                                                vacuum,               // material
                                                "logicExpHall",       // name
                                                0,                    // field manager
                                                0,                    // sensitive detector
                                                0 );                  // user limits
  fPhysExperimentalHall = new G4PVPlacement( 0,                       // rotation
                                             G4ThreeVector(),         // translation
                                             "expHall",               // name
                                             fLogicExperimentalHall,  // logical volume
                                             0,                       // mother physical volume
                                             false,                   // boolean operation
                                             0 );                     // copy number
  // Cylinder along the z-axis, with inner and outer diameter
  G4Tubs* solidTargetLayer = new G4Tubs( "solidTargetLayer",
                                         fTargetInnerRadius,          // inner radius
                                         fTargetOuterRadius,          // outer radius
                                         halfLength,                  // half cylinder length in z
                                         0.0,                         // starting phi angle in rad
                                         2.0*pi );                    // final phi angle in rad
  fLogicTargetLayer = new G4LogicalVolume( solidTargetLayer,          // solid
                                           fTargetMaterial,           // material
                                           "logicTargetLayer",        // name
                                           0,                         // field manager
                                           0,                         // sensitive detector
                                           0 );                       // user limits
  fPhysTargetLayer = new G4PVPlacement( 0,                            // rotation
                                        G4ThreeVector(),              // translation
                                        "physTargetLayer",            // name
                                        fLogicTargetLayer,            // logical volume
                                        fPhysExperimentalHall,        // mother physical volume
                                        false,                        // boolean operation
                                        0 );                          // copy number
  PrintParameters();
  G4cout << G4endl
         << "DetectorConstruction::ConstructLayer() : " << G4endl
         << "\t World (box) size: " << G4endl
         << "\t \t x : -/+ " << expHall_x << " mm ;"
         << "\t y : -/+ "    << expHall_y << " mm ;"
         << "\t z : -/+ "    << expHall_z << " mm ;" << G4endl
         << "\t Layer (cylinder) size : " << G4endl
         << "\t \t radii : " << fTargetInnerRadius << " , " << fTargetOuterRadius << " mm ;"
         << "\t \t length (along z) : " << 2.0*halfLength << " mm ;" << G4endl
         << G4endl << G4endl;
  return fPhysExperimentalHall;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial( const G4String name ) {
  fTargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial( name ); 
  if ( ! fTargetMaterial ) {
    G4cout << G4endl << G4endl
           << "WARNING: the name of the material has not been recognized!" << G4endl
           << "     ===> the default  * G4_Be *  will be used." << G4endl << G4endl;
    //***LOOKHERE*** Default material
    fTargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Be" );
  }
  if ( fLogicTargetLayer ) fLogicTargetLayer->SetMaterial( fTargetMaterial );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry() {
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  PrintParameters();
  // Update also the position of the gun
  const PrimaryGeneratorAction* pPrimaryAction = dynamic_cast< const PrimaryGeneratorAction* >
    ( G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction() );
  if ( pPrimaryAction ) pPrimaryAction->SetGunPosition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters() {
  G4cout << G4endl << G4endl
         << " ------  DetectorConstruction::PrintParameters() ------ " << G4endl
         << " Material            = " << fTargetMaterial->GetName() << G4endl
         << " Target Inner Radius = " << fTargetInnerRadius  << " mm" << G4endl
         << " Target Outer Radius = " << fTargetOuterRadius  << " mm" << G4endl
         << " B [T]               = "
         << ( fUniformMagField ? fUniformMagField->GetConstantFieldValue()/CLHEP::tesla :
              G4ThreeVector(0.0, 0.0, 0.0) ) << G4endl
         << " ------------------------------------------------------ " << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMagField( const G4double fieldValue ) {
  if ( fUniformMagField ) delete fUniformMagField;
  if ( std::abs( fieldValue ) > 0.0 ) {
    // Apply a global uniform magnetic field along the Z axis
    fUniformMagField = new G4UniformMagField( G4ThreeVector( 0.0, 0.0, fieldValue ) );
    fFieldMgr->SetDetectorField( fUniformMagField );
    fFieldMgr->CreateChordFinder( fUniformMagField );
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
