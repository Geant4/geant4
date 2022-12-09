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
           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() :  
  fMaterial( nullptr ), fExperimentalHall_log( nullptr ), fExperimentalHall_phys( nullptr ),
  fLogicLayer( nullptr ), fPhysiLayer( nullptr ),
  fLogicScoringUpDown( nullptr ), fPhysiScoringUpstream( nullptr ),
  fPhysiScoringDownstream( nullptr ),
  fLogicScoringSide( nullptr ), fPhysiScoringSide( nullptr ),
  fDetectorMessenger( nullptr ),
  fThickness( 2.0*CLHEP::m ), fDiameter( 2.0*CLHEP::m )  //***LOOKHERE*** Default values
{
  fMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Fe" );  //***LOOKHERE***
                                                                          // Default material
  fDetectorMessenger = new DetectorMessenger( this );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {
  return ConstructLayer();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructLayer() {
  // Clean old geometry, if any.
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // The target layer is a cylinder, with axis along the z-direction,
  // and positioned at the center, (0.0, 0.0, 0.0).
  // The world volume (experimental hall) is a box 20% bigger than the target layer,
  // and it is filled of "G4_Galactic" material.

  G4double expHall_x = 0.6*fDiameter;   // half dimension along x : 20% bigger than the radius
                                        //                          of the target layer 
  G4double expHall_y = 0.6*fDiameter;   // half dimension along y : 20% bigger than the radius
                                        //                          of the target layer
  G4double expHall_z = 0.6*fThickness;  // half dimension along z : 20% bigger than the half
                                        //                          thickness of the target layer

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

  // Target
  G4Tubs* solidLayer = new G4Tubs( "solidLayer",            // name
                                   0.0,                     // inner radius
                                   0.5*fDiameter,           // outer radius
                                   0.5*fThickness,          // half cylinder length in z
                                   0.0,                     // starting phi angle in rad
                                   2.0*pi );                // final phi angle in rad
  fLogicLayer = new G4LogicalVolume( solidLayer,            // solid 
                                     fMaterial,             // material
                                     "logicLayer",          // name
                                     0,                     // field manager
                                     0,                     // sensitive detector
                                     0 );                   // user limits
  fPhysiLayer = new G4PVPlacement( 0,                       // rotation
                                   G4ThreeVector(),         // translation
                                   "physiLayer",            // name
                                   fLogicLayer,             // logical volume
                                   fExperimentalHall_phys,  // mother physical volume
                                   false,                   // boolean operation
                                   0 );                     // copy number

  // Three scoring volumes: one thin layer downstream of the target ("down")
  //                        one thin layer surrounding (lateral) of the target ("side")
  //                        one thin layer upstream of the target ("up")
  G4Tubs* solidScoringUpDown = new G4Tubs( "solidScoringUpDown",    // name
                                           0.0,                     // inner radius
                                           0.5*fDiameter,           // outer radius
                                           0.5*fScoringThickness,   // half cylinder length in z
                                           0.0,                     // starting phi angle in rad
                                           2.0*pi );                // final phi angle in rad
  fLogicScoringUpDown = new G4LogicalVolume( solidScoringUpDown,    // solid 
                                             vacuum,                // material
                                             "logicScoringUpDown",  // name
                                             0,                     // field manager
                                             0,                     // sensitive detector
                                             0 );                   // user limits
  G4double zScoringUpDown = 0.5*(fThickness + fScoringThickness);
  fPhysiScoringUpstream = new G4PVPlacement( 0,                          // rotation
                                             G4ThreeVector( 0.0, 0.0, -zScoringUpDown ),
                                                                         // translation
                                             "physiScoringUpstream",     // name
                                             fLogicScoringUpDown,        // logical volume
                                             fExperimentalHall_phys,     // mother physical volume
                                             false,                      // boolean operation
                                             0 );                        // copy number  
  fPhysiScoringDownstream = new G4PVPlacement( 0,                        // rotation
                                               G4ThreeVector( 0.0, 0.0, zScoringUpDown ),
                                                                         // translation
                                               "physiScoringDownstream", // name
                                               fLogicScoringUpDown,      // logical volume
                                               fExperimentalHall_phys,   // mother physical volume
                                               false,                    // boolean operation
                                               0 );                      // copy number

  G4Tubs* solidScoringSide = new G4Tubs( "solidScoringSide",      // name
                                         0.5*fDiameter,           // inner radius
                                         0.5*fDiameter + fScoringThickness,  // outer radius
                                         0.5*fThickness,          // half cylinder length in z
                                         0.0,                     // starting phi angle in rad
                                         2.0*pi );                // final phi angle in rad
  fLogicScoringSide = new G4LogicalVolume( solidScoringSide,      // solid 
                                           vacuum,                // material
                                           "logicScoringSide",    // name
                                           0,                     // field manager
                                           0,                     // sensitive detector
                                           0 );                   // user limits
  fPhysiScoringSide = new G4PVPlacement( 0,                       // rotation
                                         G4ThreeVector( 0.0, 0.0, 0.0 ),   // translation
                                         "physiScoringSide",      // name
                                         fLogicScoringSide,       // logical volume
                                         fExperimentalHall_phys,  // mother physical volume
                                         false,                   // boolean operation
                                         0 );                     // copy number

  G4cout << G4endl
         << "DetectorConstruction::ConstructLayer() : " << G4endl
         << "\t World (box) size: " << G4endl
         << "\t \t x : -/+ " << expHall_x << " mm ;"
         << "\t y : -/+ "    << expHall_y << " mm ;"
         << "\t z : -/+ "    << expHall_z << " mm ;" << G4endl
         << "\t Target layer (cylinder) size: " << G4endl
         << "\t \t x : -/+ " << 0.5*fDiameter << " mm ;"
         << "\t y : -/+ "    << 0.5*fDiameter << " mm ;"
         << "\t z : -/+ "    << 0.5*fThickness << " mm ;" << G4endl
         << G4endl << G4endl;

  return fExperimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial( const G4String name ) {
  fMaterial = G4NistManager::Instance()->FindOrBuildMaterial( name ); 
  if ( ! fMaterial ) {
    G4cout << G4endl << G4endl
           << "WARNING: the name of the material has not been recognized!" << G4endl
           << "     ===> the default  * G4_Fe *  will be used." 
           << G4endl << G4endl;  
    fMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Fe" ); 
  }
  if ( fLogicLayer ) fLogicLayer->SetMaterial( fMaterial );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry() {
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  PrintParameters();
  // Update also the position of the gun
  const PrimaryGeneratorAction* pPrimaryAction = dynamic_cast< const PrimaryGeneratorAction* >(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction() );
  if ( pPrimaryAction ) pPrimaryAction->SetGunPosition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters() {
  G4cout << G4endl << G4endl
         << " ------  DetectorConstruction::PrintParameters() ------ " << G4endl
         << " Material         = " << fMaterial->GetName() << G4endl
         << " Thickness        = " << fThickness << " mm" << G4endl
         << " Diameter         = " << fDiameter  << " mm" << G4endl
         << " ScoringThickness = " << fScoringThickness << " mm" << G4endl 
         << " ------------------------------------------------------ "
         << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
