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
  fMaterialTracker( nullptr ), fMaterialEmCalo( nullptr ), fMaterialHadCalo( nullptr ),
  fExperimentalHall_log( nullptr ), fExperimentalHall_phys( nullptr ),
  fLogicTrackerShell( nullptr), fPhysiTrackerShell( nullptr ),
  fLogicEmCaloShell( nullptr), fPhysiEmCaloShell( nullptr),
  fLogicHadCaloShell( nullptr), fPhysiHadCaloShell( nullptr), 
  fLogicScoringTrackerShell( nullptr), fPhysiScoringTrackerShell( nullptr ),
  fLogicScoringEmCaloShell( nullptr), fPhysiScoringEmCaloShell( nullptr),
  fLogicScoringHadCaloShell( nullptr), fPhysiScoringHadCaloShell( nullptr), 
  fDetectorMessenger( nullptr ),
  fInnerRadiusTracker( 10.0*cm ), fOuterRadiusTracker(  20.0*cm ),  //***LOOKHERE*** Default radii
  fInnerRadiusEmCalo(  30.0*cm ), fOuterRadiusEmCalo(   60.0*cm ),
  fInnerRadiusHadCalo( 70.0*cm ), fOuterRadiusHadCalo( 170.0*cm )
{
  //G4cout << " BEGIN  DetectorConstruction::DetectorConstruction()" << G4endl;
  fMaterialTracker = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Si" );  //***LOOKHERE***
                                                                              // Default material
  fMaterialEmCalo  = G4NistManager::Instance()->FindOrBuildMaterial( "G4_PbWO4" );
  fMaterialHadCalo = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Fe" );
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
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructDetector() {
  //G4cout << " BEGIN  DetectorConstruction::ConstructDetector()" << G4endl;

  // Clean old geometry, if any.
  G4GeometryManager::GetInstance()->OpenGeometry();

  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Check that the radii are resonable
  G4bool isOK = true;
  if ( fInnerRadiusTracker < 0.0                                     ||
       fOuterRadiusTracker < fInnerRadiusTracker                     ||
       fInnerRadiusEmCalo  < fOuterRadiusTracker + fScoringThickness ||
       fOuterRadiusEmCalo  < fInnerRadiusEmCalo                      ||
       fInnerRadiusHadCalo < fOuterRadiusEmCalo  + fScoringThickness ||
       fOuterRadiusHadCalo < fInnerRadiusHadCalo ) {
    isOK = false;
  }
  if ( ! isOK ) {
    G4cerr << G4endl << "ERROR: the radii are inconsistent !" << G4endl
           << " InnerRadiusTracker = " << fInnerRadiusTracker << " mm" << G4endl
           << " OuterRadiusTracker = " << fOuterRadiusTracker << " mm" << G4endl
           << " InnerRadiusEmCalo  = " << fInnerRadiusEmCalo  << " mm" << G4endl
           << " OuterRadiusEmCalo  = " << fOuterRadiusEmCalo  << " mm" << G4endl
           << " InnerRadiusHadCalo = " << fInnerRadiusHadCalo << " mm" << G4endl
           << " OuterRadiusHadCalo = " << fOuterRadiusHadCalo << " mm" << G4endl
           << " ScoringThickness   = " << fScoringThickness   << " mm" << G4endl
           << G4endl;
    return nullptr;
  }
  
  // The detector consists of 3 concentric full spherical shells (G4Sphere),
  // positioned at the center, (0.0, 0.0, 0.0).
  // The world volume (experimental hall) is a box 10% bigger than the outmost
  // spherical shell.
  // and it is filled of "G4_Galactic" material.

  G4double expHall_x = 1.1*fOuterRadiusHadCalo;  // half dimension along x
  G4double expHall_y = 1.1*fOuterRadiusHadCalo;  // half dimension along y
  G4double expHall_z = 1.1*fOuterRadiusHadCalo;  // half dimension along z

  G4Material* vacuum = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Galactic" );

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

  // 1st (innermost) spherical shell: Tracker 
  G4Sphere* solidTrackerShell = new G4Sphere( "solidTrackerShell",   // name
                                              fInnerRadiusTracker,   // Inner radius
                                              fOuterRadiusTracker,   // Outer radius
                                              0.0,                   // Starting Phi angle of the
                                                                     // segment in radians
                                              2.0*CLHEP::pi,         // Delta Phi angle of the
                                                                     // segment in radians
                                              0.0,                   // Starting Theta angle of
                                                                     // the segment in radians
                                              CLHEP::pi );           // Delta Theta angle of the
                                                                     // segment in radians
  fLogicTrackerShell = new G4LogicalVolume( solidTrackerShell,       // solid 
                                            fMaterialTracker,        // material
                                            "logicTrackerShell",     // name
                                            0,                       // field manager
                                            0,                       // sensitive detector
                                            0 );                     // user limits
  fPhysiTrackerShell = new G4PVPlacement( 0,                         // rotation
                                          G4ThreeVector(),           // translation
                                          "physiTrackerShell",       // name
                                          fLogicTrackerShell,        // logical volume
                                          fExperimentalHall_phys,    // mother physical volume
                                          false,                     // boolean operation
                                          0 );                       // copy number

  // Scoring tracker shell (a thin vacuum layer, immediately outside the Tracker shell)
  G4Sphere* solidScoringTrackerShell =
    new G4Sphere( "solidScoringTrackerShell",               // name
                  fOuterRadiusTracker,                      // Inner radius
                  fOuterRadiusTracker + fScoringThickness,  // Outer radius
                  0.0,                                      // Starting Phi angle of the segment
                                                            // in radians
                  2.0*CLHEP::pi,                            // Delta Phi angle of the segment
                                                            // in radians
                  0.0,                                      // Starting Theta angle of the segment
                                                            // in radians
                  CLHEP::pi );                              // Delta Theta angle of the segment
                                                            // in radians
  fLogicScoringTrackerShell = new G4LogicalVolume( solidScoringTrackerShell,    // solid 
                                                   vacuum,                      // material
                                                   "logicScoringTrackerShell",  // name
                                                   0,                           // field manager
                                                   0,                           // sensitive
                                                                                // detector
                                                   0 );                         // user limits
  fPhysiScoringTrackerShell = new G4PVPlacement( 0,                             // rotation
                                                 G4ThreeVector(),               // translation
                                                 "physiScoringTrackerShell",    // name
                                                 fLogicScoringTrackerShell,     // logical volume
                                                 fExperimentalHall_phys,        // mother physical
                                                                                // volume
                                                 false,                         // boolean
                                                                                // operation
                                                 0 );                           // copy number
  
  // 2nd (middle) spherical shell: EM Calo 
  G4Sphere* solidEmCaloShell = new G4Sphere( "solidEmCaloShell",   // name
                                              fInnerRadiusEmCalo,  // Inner radius
                                              fOuterRadiusEmCalo,  // Outer radius
                                              0.0,                 // Starting Phi angle of the
                                                                   // segment in radians
                                              2.0*CLHEP::pi,       // Delta Phi angle of the
                                                                   // segment in radians
                                              0.0,                 // Starting Theta angle of the
                                                                   // segment in radians
                                              CLHEP::pi );         // Delta Theta angle of the
                                                                   // segment in radians
  fLogicEmCaloShell = new G4LogicalVolume( solidEmCaloShell,       // solid 
                                           fMaterialEmCalo,        // material
                                           "logicEmCaloShell",     // name
                                           0,                      // field manager
                                           0,                      // sensitive detector
                                           0 );                    // user limits
  fPhysiEmCaloShell = new G4PVPlacement( 0,                        // rotation
                                         G4ThreeVector(),          // translation
                                         "physiEmCaloShell",       // name
                                         fLogicEmCaloShell,        // logical volume
                                         fExperimentalHall_phys,   // mother physical volume
                                         false,                    // boolean operation
                                         0 );                      // copy number

  // Scoring EmCalo shell (a thin vacuum layer, immediately outside the EmCalo shell)
  G4Sphere* solidScoringEmCaloShell =
    new G4Sphere( "solidScoringEmCaloShell",               // name
                  fOuterRadiusEmCalo,                      // Inner radius
                  fOuterRadiusEmCalo + fScoringThickness,  // Outer radius
                  0.0,                                     // Starting Phi angle of the segment
                                                           // in radians
                  2.0*CLHEP::pi,                           // Delta Phi angle of the segment
                                                           // in radians
                  0.0,                                     // Starting Theta angle of the
                                                           // segment in radians
                  CLHEP::pi );                             // Delta Theta angle of the segment
                                                           // in radians
  fLogicScoringEmCaloShell = new G4LogicalVolume( solidScoringEmCaloShell,    // solid 
                                                  vacuum,                     // material
                                                  "logicScoringEmCaloShell",  // name
                                                  0,                          // field manager
                                                  0,                          // sensitive
                                                                              // detector
                                                  0 );                        // user limits
  fPhysiScoringEmCaloShell = new G4PVPlacement( 0,                            // rotation
                                                G4ThreeVector(),              // translation
                                                "physiScoringEmCaloShell",    // name
                                                fLogicScoringEmCaloShell,     // logical volume
                                                fExperimentalHall_phys,       // mother physical
                                                                              // volume
                                                false,                        // boolean operation
                                                0 );                          // copy number
  
  // 3rd (outmost) spherical shell: HAD Calo 
  G4Sphere* solidHadCaloShell = new G4Sphere( "solidHadCaloShell",  // name
                                              fInnerRadiusHadCalo,  // Inner radius
                                              fOuterRadiusHadCalo,  // Outer radius
                                              0.0,                  // Starting Phi angle of the
                                                                    // segment in radians
                                              2.0*CLHEP::pi,        // Delta Phi angle of the
                                                                    // segment in radians
                                              0.0,                  // Starting Theta angle of
                                                                    // the segment in radians
                                              CLHEP::pi );          // Delta Theta angle of the
                                                                    // segment in radians
  fLogicHadCaloShell = new G4LogicalVolume( solidHadCaloShell,      // solid 
                                            fMaterialHadCalo,       // material
                                            "logicHadCaloShell",    // name
                                            0,                      // field manager
                                            0,                      // sensitive detector
                                            0 );                    // user limits
  fPhysiHadCaloShell = new G4PVPlacement( 0,                        // rotation
                                          G4ThreeVector(),          // translation
                                          "physiHadCaloShell",      // name
                                          fLogicHadCaloShell,       // logical volume
                                          fExperimentalHall_phys,   // mother physical volume
                                          false,                    // boolean operation
                                          0 );                      // copy number

  // Scoring HadCalo shell (a thin vacuum layer, immediately outside the HadCalo shell)
  G4Sphere* solidScoringHadCaloShell =
    new G4Sphere( "solidScoringHadCaloShell",               // name
                  fOuterRadiusHadCalo,                      // Inner radius
                  fOuterRadiusHadCalo + fScoringThickness,  // Outer radius
                  0.0,                                      // Starting Phi angle of the
                                                            // segment in radians
                  2.0*CLHEP::pi,                            // Delta Phi angle of the segment
                                                            // in radians
                  0.0,                                      // Starting Theta angle of the
                                                            // segment in radians
                  CLHEP::pi );                              // Delta Theta angle of the segment
                                                            // in radians
  fLogicScoringHadCaloShell = new G4LogicalVolume( solidScoringHadCaloShell,    // solid 
                                                   vacuum,                      // material
                                                   "logicScoringHadCaloShell",  // name
                                                   0,                           // field manager
                                                   0,                           // sensitive
                                                                                // detector
                                                   0 );                         // user limits
  fPhysiScoringHadCaloShell = new G4PVPlacement( 0,                             // rotation
                                                 G4ThreeVector(),               // translation
                                                 "physiScoringHadCaloShell",    // name
                                                 fLogicScoringHadCaloShell,     // logical volume
                                                 fExperimentalHall_phys,        // mother physical
                                                                                // volume
                                                 false,                         // boolean
                                                                                // operation
                                                 0 );                           // copy number
  
  G4cout << "DetectorConstruction::ConstructSphere() : " << G4endl
         << "\t World (box) size: " << G4endl
         << "\t \t x : -/+ " << expHall_x << " mm ;"
         << "\t y : -/+ "    << expHall_y << " mm ;"
         << "\t z : -/+ "    << expHall_z << " mm ;" << G4endl
         << G4endl;
  PrintParameters();
  //G4cout << " END  DetectorConstruction::ConstructDetector()

  return fExperimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterialTracker( const G4String name ) {
  fMaterialTracker = G4NistManager::Instance()->FindOrBuildMaterial( name ); 
  if ( ! fMaterialTracker ) {
    G4cout << G4endl << G4endl
           << "WARNING: the name of the material has not been recognized!" << G4endl
           << "     ===> the default  * G4_Si *  will be used." 
           << G4endl << G4endl;  
    fMaterialTracker = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Si" ); 
  }
  if ( fLogicTrackerShell ) fLogicTrackerShell->SetMaterial( fMaterialTracker );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterialEmCalo( const G4String name ) {
  fMaterialEmCalo = G4NistManager::Instance()->FindOrBuildMaterial( name ); 
  if ( ! fMaterialEmCalo ) {
    G4cout << G4endl << G4endl
           << "WARNING: the name of the material has not been recognized!" << G4endl
           << "     ===> the default  * G4_Pb *  will be used." 
           << G4endl << G4endl;  
    fMaterialEmCalo = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Pb" ); 
  }
  if ( fLogicEmCaloShell ) fLogicEmCaloShell->SetMaterial( fMaterialEmCalo );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterialHadCalo( const G4String name ) {
  fMaterialHadCalo = G4NistManager::Instance()->FindOrBuildMaterial( name ); 
  if ( fMaterialHadCalo == nullptr ) {
    G4cout << G4endl << G4endl
           << "WARNING: the name of the material has not been recognized!" << G4endl
           << "     ===> the default  * G4_Fe *  will be used." 
           << G4endl << G4endl;  
    fMaterialHadCalo = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Fe" ); 
  }
  if ( fLogicHadCaloShell ) fLogicHadCaloShell->SetMaterial( fMaterialHadCalo );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry() {
  //G4cout << " BEGIN  DetectorConstruction::UpdateGeometry" << G4endl;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  PrintParameters();
  // Update also the position of the gun
  const PrimaryGeneratorAction* pPrimaryAction = 
    dynamic_cast< const PrimaryGeneratorAction* >
      ( G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction() );
  if ( pPrimaryAction ) pPrimaryAction->SetGunPosition();
  //G4cout << " END  DetectorConstruction::UpdateGeometry" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters() {
  G4cout << G4endl
         << " ------  DetectorConstruction::PrintParameters() ------ " << G4endl
         << " MaterialTracker = " << fMaterialTracker->GetName() << G4endl
         << " MaterialEmCalo  = " << fMaterialEmCalo->GetName() << G4endl
         << " MaterialHadCalo = " << fMaterialHadCalo->GetName() << G4endl
         << " InnerRadiusTracker = " << fInnerRadiusTracker << " mm" << G4endl
         << " OuterRadiusTracker = " << fOuterRadiusTracker << " mm" << G4endl
         << " InnerRadiusEmCalo  = " << fInnerRadiusEmCalo  << " mm" << G4endl
         << " OuterRadiusEmCalo  = " << fOuterRadiusEmCalo  << " mm" << G4endl
         << " InnerRadiusHadCalo = " << fInnerRadiusHadCalo << " mm" << G4endl
         << " OuterRadiusHadCalo = " << fOuterRadiusHadCalo << " mm" << G4endl
         << " ScoringThickness   = " << fScoringThickness   << " mm" << G4endl
         << " -------------------------------------------------------- " << G4endl
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
