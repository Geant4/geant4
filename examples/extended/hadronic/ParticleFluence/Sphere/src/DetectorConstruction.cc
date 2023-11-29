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
#include "G4Orb.hh"
#include "G4Sphere.hh"
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
           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() :
  fMaterial( nullptr ),
  fExperimentalHall_log( nullptr ), fExperimentalHall_phys( nullptr ),
  fLogicSphere( nullptr ), fPhysiSphere( nullptr ),
  fLogicScoringShell( nullptr ), fPhysiScoringShell( nullptr ),
  fDetectorMessenger( nullptr ),
  fRadius( 1.0*CLHEP::m )  //***LOOKHERE*** Default values
{
  //G4cout << " BEGIN  DetectorConstruction::DetectorConstruction()" << G4endl;
  fMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Fe" );  //***LOOKHERE***
                                                                          // Default material
  fDetectorMessenger = new DetectorMessenger( this );
  //G4cout << " END  DetectorConstruction::DetectorConstruction()" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {
  //G4cout << " BEGIN  DetectorConstruction::Construct()" << G4endl;
  return ConstructSphere();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructSphere() {
  //G4cout << " BEGIN  DetectorConstruction::ConstructSphere()" << G4endl;

  // Clean old geometry, if any.
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // The target is a full solid sphere (G4Orb), positioned at the center, (0.0, 0.0, 0.0).
  // The world volume (experimental hall) is a box 10% bigger than the sphere
  // and it is filled of "G4_Galactic" material.

  G4double expHall_x = 1.1*fRadius;  // half dimension along x
  G4double expHall_y = 1.1*fRadius;  // half dimension along y
  G4double expHall_z = 1.1*fRadius;  // half dimension along z

  G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Galactic" );

  // Experimental hall
  G4Box* experimentalHall_box = new G4Box( "expHall_box", expHall_x, expHall_y, expHall_z );
  fExperimentalHall_log = new G4LogicalVolume( experimentalHall_box, // solid 
                                               vacuum,               // material
                                               "expHall_log",        // name
                                               0,                    // field manager
                                               0,                    // sensitive detector
                                               0 );                  // user limits
  fExperimentalHall_phys = new G4PVPlacement( 0,                     // rotation
                                              G4ThreeVector(),       // translation
                                              "expHall",             // name
                                              fExperimentalHall_log, // logical volume
                                              0,                     // mother physical volume
                                              false,                 // boolean operation
                                              0 );                   // copy number

  // Target sphere
  G4Orb* solidSphere = new G4Orb( "solidSphere",            // name
                                  fRadius );                // outer radius
  fLogicSphere = new G4LogicalVolume( solidSphere,          // solid
                                      fMaterial,            // material
                                      "logicSphere",        // name
                                      0,                    // field manager
                                      0,                    // sensitive detector
                                      0 );                  // user limits
  fPhysiSphere = new G4PVPlacement( 0,                      // rotation
                                    G4ThreeVector(),        // translation
                                    "physiSphere",          // name
                                    fLogicSphere,           // logical volume
                                    fExperimentalHall_phys, // mother physical volume
                                    false,                  // boolean operation
                                    0 );                    // copy number

  // Scoring shell (a thin vacuum layer, immediately outside the target sphere)
  G4Sphere* solidScoringShell = new G4Sphere( "solidScoringShell",   // name
                                              fRadius,               // Inner radius (the radius
                                                                     // of the target sphere)
                                              fRadius + fScoringThickness,  // Outer radius
                                              0.0,                   // Starting Phi angle of the
                                                                     // segment in radians
                                              2.0*CLHEP::pi,         // Delta Phi angle of the
                                                                     // segment in radians
                                              0.0,                   // Starting Theta angle of
                                                                     // the segment in radians
                                              CLHEP::pi );           // Delta Theta angle of the
                                                                     // segment in radians
  fLogicScoringShell = new G4LogicalVolume( solidScoringShell,       // solid
                                            vacuum,                  // material
                                            "logicScoringShell",     // name
                                            0,                       // field manager
                                            0,                       // sensitive detector
                                            0 );                     // user limits
  fPhysiScoringShell = new G4PVPlacement( 0,                         // rotation
                                          G4ThreeVector(),           // translation
                                          "physiScoringShell",       // name
                                          fLogicScoringShell,        // logical volume
                                          fExperimentalHall_phys,    // mother physical volume
                                          false,                     // boolean operation
                                          0 );                       // copy number

  G4cout << G4endl
         << "DetectorConstruction::ConstructSphere() : " << G4endl
         << "\t World (box) size: " << G4endl
         << "\t \t x : -/+ " << expHall_x << " mm ;"
         << "\t y : -/+ "    << expHall_y << " mm ;"
         << "\t z : -/+ "    << expHall_z << " mm ;" << G4endl
         << G4endl << G4endl;

  //G4cout << " END  DetectorConstruction::ConstructSphere()

  return fExperimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial( const G4String name ) {
  fMaterial = G4NistManager::Instance()->FindOrBuildMaterial( name );
  if ( fMaterial == nullptr ) {
    G4cout << G4endl << G4endl
           << "WARNING: the name of the material has not been recognized!" << G4endl
           << "     ===> the default  * G4_Fe *  will be used."
           << G4endl << G4endl;
    fMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Fe" );
  }
  if ( fLogicSphere ) fLogicSphere->SetMaterial( fMaterial );
  //G4cout << " Absorber Material = " << logicSphere->GetMaterial()->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry() {
  //G4cout << " BEGIN  DetectorConstruction::UpdateGeometry" << G4endl;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  PrintParameters();
  // Update also the position of the gun
  const PrimaryGeneratorAction* pPrimaryAction = dynamic_cast< const PrimaryGeneratorAction* >
    ( G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction() );
  if ( pPrimaryAction ) pPrimaryAction->SetGunPosition();
  //G4cout << " END  DetectorConstruction::UpdateGeometry" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters() {
  G4cout << G4endl << G4endl
         << " ------  DetectorConstruction::PrintParameters() ------ " << G4endl
         << " Material         = " << fMaterial->GetName() << G4endl
         << " Radius           = " << fRadius << " mm" << G4endl
         << " ScoringThickness = " << fScoringThickness << " mm" << G4endl
         << " -------------------------------------------------------- "
         << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
