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
// $Id: B5DetectorConstruction.cc 103284 2017-03-24 08:27:19Z gcosmo $
//
/// \file B5DetectorConstruction.cc
/// \brief Implementation of the B5DetectorConstruction class

#include "B5DetectorConstruction.hh"
#include "B5MagneticField.hh"
#include "B5CellParameterisation.hh"
#include "B5HodoscopeSD.hh"
#include "B5DriftChamberSD.hh"
#include "B5EmCalorimeterSD.hh"
#include "B5HadCalorimeterSD.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal B5MagneticField* B5DetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* B5DetectorConstruction::fFieldMgr = 0;
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DetectorConstruction::B5DetectorConstruction()
: G4VUserDetectorConstruction(), 
  fMessenger(nullptr),
  fHodoscope1Logical(nullptr), fHodoscope2Logical(nullptr),
  fWirePlane1Logical(nullptr), fWirePlane2Logical(nullptr),
  fCellLogical(nullptr), fHadCalScintiLogical(nullptr),
  fMagneticLogical(nullptr),
  fVisAttributes(),
  fArmAngle(30.*deg), fArmRotation(nullptr), fSecondArmPhys(nullptr)

{
  fArmRotation = new G4RotationMatrix();
  fArmRotation->rotateY(fArmAngle);
  
  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DetectorConstruction::~B5DetectorConstruction()
{
  delete fArmRotation;
  delete fMessenger;
  
  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B5DetectorConstruction::Construct()
{
  // Construct materials
  ConstructMaterials();
  auto air = G4Material::GetMaterial("G4_AIR");
  //auto argonGas = G4Material::GetMaterial("B5_Ar");
  auto argonGas = G4Material::GetMaterial("G4_Ar");
  auto scintillator = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto csI = G4Material::GetMaterial("G4_CESIUM_IODIDE");
  auto lead = G4Material::GetMaterial("G4_Pb");
  
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // geometries --------------------------------------------------------------
  // experimental hall (world volume)
  auto worldSolid 
    = new G4Box("worldBox",10.*m,3.*m,10.*m);
  auto worldLogical
    = new G4LogicalVolume(worldSolid,air,"worldLogical");
  auto worldPhysical
    = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
                        false,0,checkOverlaps);
  
  // Tube with Local Magnetic field
  
  auto magneticSolid 
    = new G4Tubs("magneticTubs",0.,1.*m,1.*m,0.,360.*deg);

  fMagneticLogical
    = new G4LogicalVolume(magneticSolid, air, "magneticLogical");

  // placement of Tube
  
  G4RotationMatrix* fieldRot = new G4RotationMatrix();
  fieldRot->rotateX(90.*deg);
  new G4PVPlacement(fieldRot,G4ThreeVector(),fMagneticLogical,
                    "magneticPhysical",worldLogical,
                    false,0,checkOverlaps);
  
  // set step limit in tube with magnetic field  
  G4UserLimits* userLimits = new G4UserLimits(1*m);
  fMagneticLogical->SetUserLimits(userLimits);
  
  // first arm
  auto firstArmSolid 
    = new G4Box("firstArmBox",1.5*m,1.*m,3.*m);
  auto firstArmLogical
    = new G4LogicalVolume(firstArmSolid,air,"firstArmLogical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,-5.*m),firstArmLogical,
                    "firstArmPhysical",worldLogical,
                    false,0,checkOverlaps);
  
  // second arm
  auto secondArmSolid 
    = new G4Box("secondArmBox",2.*m,2.*m,3.5*m);
  auto secondArmLogical
    = new G4LogicalVolume(secondArmSolid,air,"secondArmLogical");
  auto x = -5.*m * std::sin(fArmAngle);
  auto z = 5.*m * std::cos(fArmAngle);
  fSecondArmPhys
    = new G4PVPlacement(fArmRotation,G4ThreeVector(x,0.,z),secondArmLogical,
                        "fSecondArmPhys",worldLogical,
                        false,0,checkOverlaps);
  
  // hodoscopes in first arm
  auto hodoscope1Solid 
    = new G4Box("hodoscope1Box",5.*cm,20.*cm,0.5*cm);
  fHodoscope1Logical
    = new G4LogicalVolume(hodoscope1Solid,scintillator,"hodoscope1Logical");

  for (auto i=0;i<kNofHodoscopes1;i++) {
      G4double x1 = (i-kNofHodoscopes1/2)*10.*cm;
      new G4PVPlacement(0,G4ThreeVector(x1,0.,-1.5*m),fHodoscope1Logical,
                        "hodoscope1Physical",firstArmLogical,
                        false,i,checkOverlaps);
  }
  
  // drift chambers in first arm
  auto chamber1Solid 
    = new G4Box("chamber1Box",1.*m,30.*cm,1.*cm);
  auto chamber1Logical
    = new G4LogicalVolume(chamber1Solid,argonGas,"chamber1Logical");

  for (auto i=0;i<kNofChambers;i++) {
    G4double z1 = (i-kNofChambers/2)*0.5*m;
    new G4PVPlacement(0,G4ThreeVector(0.,0.,z1),chamber1Logical,
                      "chamber1Physical",firstArmLogical,
                      false,i,checkOverlaps);
  }
  
  // "virtual" wire plane
  auto wirePlane1Solid 
    = new G4Box("wirePlane1Box",1.*m,30.*cm,0.1*mm);
  fWirePlane1Logical
    = new G4LogicalVolume(wirePlane1Solid,argonGas,"wirePlane1Logical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fWirePlane1Logical,
                    "wirePlane1Physical",chamber1Logical,
                    false,0,checkOverlaps);
  
  // hodoscopes in second arm
  auto hodoscope2Solid 
    = new G4Box("hodoscope2Box",5.*cm,20.*cm,0.5*cm);
  fHodoscope2Logical
    = new G4LogicalVolume(hodoscope2Solid,scintillator,"hodoscope2Logical");

  for (auto i=0;i<kNofHodoscopes2;i++) {
      G4double x2 = (i-kNofHodoscopes2/2)*10.*cm;
      new G4PVPlacement(0,G4ThreeVector(x2,0.,0.),fHodoscope2Logical,
                        "hodoscope2Physical",secondArmLogical,
                        false,i,checkOverlaps);
  }
  
  // drift chambers in second arm
  auto chamber2Solid 
    = new G4Box("chamber2Box",1.5*m,30.*cm,1.*cm);
  auto chamber2Logical
    = new G4LogicalVolume(chamber2Solid,argonGas,"chamber2Logical");
  
  for (auto i=0;i<kNofChambers;i++) {
    G4double z2 = (i-kNofChambers/2)*0.5*m - 1.5*m;
    new G4PVPlacement(0,G4ThreeVector(0.,0.,z2),chamber2Logical,
                      "chamber2Physical",secondArmLogical,
                      false,i,checkOverlaps);
  }
  
  // "virtual" wire plane
  auto wirePlane2Solid 
    = new G4Box("wirePlane2Box",1.5*m,30.*cm,0.1*mm);
  fWirePlane2Logical
    = new G4LogicalVolume(wirePlane2Solid,argonGas,"wirePlane2Logical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fWirePlane2Logical,
                    "wirePlane2Physical",chamber2Logical,
                    false,0,checkOverlaps);
  
  // CsI calorimeter
  auto emCalorimeterSolid 
    = new G4Box("EMcalorimeterBox",1.5*m,30.*cm,15.*cm);
  auto emCalorimeterLogical
    = new G4LogicalVolume(emCalorimeterSolid,csI,"EMcalorimeterLogical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,2.*m),emCalorimeterLogical,
                    "EMcalorimeterPhysical",secondArmLogical,
                    false,0,checkOverlaps);
  
  // EMcalorimeter cells
  auto cellSolid 
    = new G4Box("cellBox",7.5*cm,7.5*cm,15.*cm);
  fCellLogical
    = new G4LogicalVolume(cellSolid,csI,"cellLogical");
  G4VPVParameterisation* cellParam = new B5CellParameterisation();
  new G4PVParameterised("cellPhysical",fCellLogical,emCalorimeterLogical,
                        kXAxis,kNofEmCells,cellParam);
  
  // hadron calorimeter
  auto hadCalorimeterSolid
    = new G4Box("HadCalorimeterBox",1.5*m,30.*cm,50.*cm);
  auto hadCalorimeterLogical
    = new G4LogicalVolume(hadCalorimeterSolid,lead,"HadCalorimeterLogical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,3.*m),hadCalorimeterLogical,
                    "HadCalorimeterPhysical",secondArmLogical,
                    false,0,checkOverlaps);
  
  // hadron calorimeter column
  auto HadCalColumnSolid
    = new G4Box("HadCalColumnBox",15.*cm,30.*cm,50.*cm);
  auto HadCalColumnLogical
    = new G4LogicalVolume(HadCalColumnSolid,lead,"HadCalColumnLogical");
  new G4PVReplica("HadCalColumnPhysical",HadCalColumnLogical,
                  hadCalorimeterLogical,kXAxis,kNofHadColumns,30.*cm);
  
  // hadron calorimeter cell
  auto HadCalCellSolid
    = new G4Box("HadCalCellBox",15.*cm,15.*cm,50.*cm);
  auto HadCalCellLogical
    = new G4LogicalVolume(HadCalCellSolid,lead,"HadCalCellLogical");
  new G4PVReplica("HadCalCellPhysical",HadCalCellLogical,
                  HadCalColumnLogical,kYAxis,kNofHadRows,30.*cm);
  
  // hadron calorimeter layers
  auto HadCalLayerSolid
    = new G4Box("HadCalLayerBox",15.*cm,15.*cm,2.5*cm);
  auto HadCalLayerLogical
    = new G4LogicalVolume(HadCalLayerSolid,lead,"HadCalLayerLogical");
  new G4PVReplica("HadCalLayerPhysical",HadCalLayerLogical,
                  HadCalCellLogical,kZAxis,kNofHadCells,5.*cm);
  
  // scintillator plates
  auto HadCalScintiSolid
    = new G4Box("HadCalScintiBox",15.*cm,15.*cm,0.5*cm);
  fHadCalScintiLogical
    = new G4LogicalVolume(HadCalScintiSolid,scintillator,
                          "HadCalScintiLogical");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,2.*cm),fHadCalScintiLogical,
                    "HadCalScintiPhysical",HadCalLayerLogical,
                    false,0,checkOverlaps);
  
  // visualization attributes ------------------------------------------------
  
  auto visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAttributes->SetVisibility(false);
  worldLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.9));   // LightGray
  fMagneticLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAttributes->SetVisibility(false);
  firstArmLogical->SetVisAttributes(visAttributes);
  secondArmLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.8888,0.0,0.0));
  fHodoscope1Logical->SetVisAttributes(visAttributes);
  fHodoscope2Logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  chamber1Logical->SetVisAttributes(visAttributes);
  chamber2Logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.0,0.8888,0.0));
  visAttributes->SetVisibility(false);
  fWirePlane1Logical->SetVisAttributes(visAttributes);
  fWirePlane2Logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.8888,0.8888,0.0));
  visAttributes->SetVisibility(false);
  emCalorimeterLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.0));
  fCellLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
  hadCalorimeterLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
  visAttributes->SetVisibility(false);
  HadCalColumnLogical->SetVisAttributes(visAttributes);
  HadCalCellLogical->SetVisAttributes(visAttributes);
  HadCalLayerLogical->SetVisAttributes(visAttributes);
  fHadCalScintiLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  // return the world physical volume ----------------------------------------
  
  return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::ConstructSDandField()
{
  // sensitive detectors -----------------------------------------------------
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;
  
  auto hodoscope1 = new B5HodoscopeSD(SDname="/hodoscope1");
  sdManager->AddNewDetector(hodoscope1);
  fHodoscope1Logical->SetSensitiveDetector(hodoscope1);

  auto hodoscope2 = new B5HodoscopeSD(SDname="/hodoscope2");
  sdManager->AddNewDetector(hodoscope2);
  fHodoscope2Logical->SetSensitiveDetector(hodoscope2);
  
  auto chamber1 = new B5DriftChamberSD(SDname="/chamber1");
  sdManager->AddNewDetector(chamber1);
  fWirePlane1Logical->SetSensitiveDetector(chamber1);

  auto chamber2 = new B5DriftChamberSD(SDname="/chamber2");
  sdManager->AddNewDetector(chamber2);
  fWirePlane2Logical->SetSensitiveDetector(chamber2);
  
  auto emCalorimeter = new B5EmCalorimeterSD(SDname="/EMcalorimeter");
  sdManager->AddNewDetector(emCalorimeter);
  fCellLogical->SetSensitiveDetector(emCalorimeter);
  
  auto hadCalorimeter = new B5HadCalorimeterSD(SDname="/HadCalorimeter");
  sdManager->AddNewDetector(hadCalorimeter);
  fHadCalScintiLogical->SetSensitiveDetector(hadCalorimeter);

  // magnetic field ----------------------------------------------------------
  fMagneticField = new B5MagneticField();
  fFieldMgr = new G4FieldManager();
  fFieldMgr->SetDetectorField(fMagneticField);
  fFieldMgr->CreateChordFinder(fMagneticField);
  G4bool forceToAllDaughters = true;
  fMagneticLogical->SetFieldManager(fFieldMgr, forceToAllDaughters);
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::ConstructMaterials()
{
  auto nistManager = G4NistManager::Instance();

  // Air 
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  // Argon gas
  nistManager->FindOrBuildMaterial("G4_Ar");
  // With a density different from the one defined in NIST
  // G4double density = 1.782e-03*g/cm3; 
  // nistManager->BuildMaterialWithNewDensity("B5_Ar","G4_Ar",density);
  // !! cases segmentation fault

  // Scintillator
  // (PolyVinylToluene, C_9H_10)
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
  // CsI
  nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  
  // Lead
  nistManager->FindOrBuildMaterial("G4_Pb");
  
  // Vacuum "Galactic"
  // nistManager->FindOrBuildMaterial("G4_Galactic");

  // Vacuum "Air with low density"
  // auto air = G4Material::GetMaterial("G4_AIR");
  // G4double density = 1.0e-5*air->GetDensity();
  // nistManager
  //   ->BuildMaterialWithNewDensity("Air_lowDensity", "G4_AIR", density);

  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::SetArmAngle(G4double val)
{
  if (!fSecondArmPhys) {
      G4cerr << "Detector has not yet been constructed." << G4endl;
      return;
  }
  
  fArmAngle = val;
  *fArmRotation = G4RotationMatrix();  // make it unit vector
  fArmRotation->rotateY(fArmAngle);
  auto x = -5.*m * std::sin(fArmAngle);
  auto z = 5.*m * std::cos(fArmAngle);
  fSecondArmPhys->SetTranslation(G4ThreeVector(x,0.,z));
  
  // tell G4RunManager that we change the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::DefineCommands()
{
  // Define /B5/detector command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, 
                                      "/B5/detector/", 
                                      "Detector control");

  // armAngle command
  auto& armAngleCmd
    = fMessenger->DeclareMethodWithUnit("armAngle","deg",
                                &B5DetectorConstruction::SetArmAngle, 
                                "Set rotation angle of the second arm.");
  armAngleCmd.SetParameterName("angle", true);
  armAngleCmd.SetRange("angle>=0. && angle<180.");
  armAngleCmd.SetDefaultValue("30.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
