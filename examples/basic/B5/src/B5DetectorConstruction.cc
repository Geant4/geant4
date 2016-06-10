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
// $Id: B5DetectorConstruction.cc 77656 2013-11-27 08:52:57Z gcosmo $
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
#include "G4AutoDelete.hh"

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
  fMessenger(0),
  fHodoscope1Logical(0), fHodoscope2Logical(0),
  fWirePlane1Logical(0), fWirePlane2Logical(0),
  fCellLogical(0), fHadCalScintiLogical(0),
  fMagneticLogical(0),
  fVisAttributes(),
  fArmAngle(30.*deg), fArmRotation(0), fSecondArmPhys(0)

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
    
    for (G4int i=0; i<G4int(fVisAttributes.size()); ++i) 
    {
      delete fVisAttributes[i];
    }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B5DetectorConstruction::Construct()
{
    // Construct materials
    ConstructMaterials();
    G4Material* air = G4Material::GetMaterial("G4_AIR");
    //G4Material* argonGas = G4Material::GetMaterial("B5_Ar");
    G4Material* argonGas = G4Material::GetMaterial("G4_Ar");
    G4Material* scintillator 
      = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4Material* csI = G4Material::GetMaterial("G4_CESIUM_IODIDE");
    G4Material* lead = G4Material::GetMaterial("G4_Pb");
    
    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = true;

    // geometries --------------------------------------------------------------
    // experimental hall (world volume)
    G4VSolid* worldSolid 
      = new G4Box("worldBox",10.*m,3.*m,10.*m);
    G4LogicalVolume* worldLogical
      = new G4LogicalVolume(worldSolid,air,"worldLogical");
    G4VPhysicalVolume* worldPhysical
      = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
                          false,0,checkOverlaps);
    
    // Tube with Local Magnetic field
    
    G4VSolid* magneticSolid 
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
    G4VSolid* firstArmSolid 
      = new G4Box("firstArmBox",1.5*m,1.*m,3.*m);
    G4LogicalVolume* firstArmLogical
      = new G4LogicalVolume(firstArmSolid,air,"firstArmLogical");
    new G4PVPlacement(0,G4ThreeVector(0.,0.,-5.*m),firstArmLogical,
                      "firstArmPhysical",worldLogical,
                      false,0,checkOverlaps);
    
    // second arm
    G4VSolid* secondArmSolid 
      = new G4Box("secondArmBox",2.*m,2.*m,3.5*m);
    G4LogicalVolume* secondArmLogical
      = new G4LogicalVolume(secondArmSolid,air,"secondArmLogical");
    G4double x = -5.*m * std::sin(fArmAngle);
    G4double z = 5.*m * std::cos(fArmAngle);
    fSecondArmPhys
      = new G4PVPlacement(fArmRotation,G4ThreeVector(x,0.,z),secondArmLogical,
                          "fSecondArmPhys",worldLogical,
                          false,0,checkOverlaps);
    
    // hodoscopes in first arm
    G4VSolid* hodoscope1Solid 
      = new G4Box("hodoscope1Box",5.*cm,20.*cm,0.5*cm);
    fHodoscope1Logical
      = new G4LogicalVolume(hodoscope1Solid,scintillator,"hodoscope1Logical");
    for (G4int i=0;i<15;i++)
    {
        G4double x1 = (i-7)*10.*cm;
        new G4PVPlacement(0,G4ThreeVector(x1,0.,-1.5*m),fHodoscope1Logical,
                          "hodoscope1Physical",firstArmLogical,
                          false,i,checkOverlaps);
    }
    
    // drift chambers in first arm
    G4VSolid* chamber1Solid 
      = new G4Box("chamber1Box",1.*m,30.*cm,1.*cm);
    G4LogicalVolume* chamber1Logical
      = new G4LogicalVolume(chamber1Solid,argonGas,"chamber1Logical");
    for (G4int i=0;i<5;i++)
    {
        G4double z1 = (i-2)*0.5*m;
        new G4PVPlacement(0,G4ThreeVector(0.,0.,z1),chamber1Logical,
                          "chamber1Physical",firstArmLogical,
                          false,i,checkOverlaps);
    }
    
    // "virtual" wire plane
    G4VSolid* wirePlane1Solid 
      = new G4Box("wirePlane1Box",1.*m,30.*cm,0.1*mm);
    fWirePlane1Logical
      = new G4LogicalVolume(wirePlane1Solid,argonGas,"wirePlane1Logical");
    new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fWirePlane1Logical,
                      "wirePlane1Physical",chamber1Logical,
                      false,0,checkOverlaps);
    
    // hodoscopes in second arm
    G4VSolid* hodoscope2Solid 
      = new G4Box("hodoscope2Box",5.*cm,20.*cm,0.5*cm);
    fHodoscope2Logical
      = new G4LogicalVolume(hodoscope2Solid,scintillator,"hodoscope2Logical");
    for (G4int i=0;i<25;i++)
    {
        G4double x2 = (i-12)*10.*cm;
        new G4PVPlacement(0,G4ThreeVector(x2,0.,0.),fHodoscope2Logical,
                          "hodoscope2Physical",secondArmLogical,
                          false,i,checkOverlaps);
    }
    
    // drift chambers in second arm
    G4VSolid* chamber2Solid 
      = new G4Box("chamber2Box",1.5*m,30.*cm,1.*cm);
    G4LogicalVolume* chamber2Logical
      = new G4LogicalVolume(chamber2Solid,argonGas,"chamber2Logical");
    for (G4int i=0;i<5;i++)
    {
        G4double z2 = (i-2)*0.5*m - 1.5*m;
        new G4PVPlacement(0,G4ThreeVector(0.,0.,z2),chamber2Logical,
                          "chamber2Physical",secondArmLogical,
                          false,i,checkOverlaps);
    }
    
    // "virtual" wire plane
    G4VSolid* wirePlane2Solid 
      = new G4Box("wirePlane2Box",1.5*m,30.*cm,0.1*mm);
    fWirePlane2Logical
      = new G4LogicalVolume(wirePlane2Solid,argonGas,"wirePlane2Logical");
    new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fWirePlane2Logical,
                      "wirePlane2Physical",chamber2Logical,
                      false,0,checkOverlaps);
    
    // CsI calorimeter
    G4VSolid* emCalorimeterSolid 
      = new G4Box("EMcalorimeterBox",1.5*m,30.*cm,15.*cm);
    G4LogicalVolume* emCalorimeterLogical
      = new G4LogicalVolume(emCalorimeterSolid,csI,"EMcalorimeterLogical");
    new G4PVPlacement(0,G4ThreeVector(0.,0.,2.*m),emCalorimeterLogical,
                      "EMcalorimeterPhysical",secondArmLogical,
                      false,0,checkOverlaps);
    
    // EMcalorimeter cells
    G4VSolid* cellSolid 
      = new G4Box("cellBox",7.5*cm,7.5*cm,15.*cm);
    fCellLogical
      = new G4LogicalVolume(cellSolid,csI,"cellLogical");
    G4VPVParameterisation* cellParam = new B5CellParameterisation();
    new G4PVParameterised("cellPhysical",fCellLogical,emCalorimeterLogical,
                          kXAxis,80,cellParam);
    
    // hadron calorimeter
    G4VSolid* hadCalorimeterSolid
      = new G4Box("HadCalorimeterBox",1.5*m,30.*cm,50.*cm);
    G4LogicalVolume* hadCalorimeterLogical
      = new G4LogicalVolume(hadCalorimeterSolid,lead,"HadCalorimeterLogical");
    new G4PVPlacement(0,G4ThreeVector(0.,0.,3.*m),hadCalorimeterLogical,
                      "HadCalorimeterPhysical",secondArmLogical,
                      false,0,checkOverlaps);
    
    // hadron calorimeter column
    G4VSolid* HadCalColumnSolid
      = new G4Box("HadCalColumnBox",15.*cm,30.*cm,50.*cm);
    G4LogicalVolume* HadCalColumnLogical
      = new G4LogicalVolume(HadCalColumnSolid,lead,"HadCalColumnLogical");
    new G4PVReplica("HadCalColumnPhysical",HadCalColumnLogical,
                    hadCalorimeterLogical,kXAxis,10,30.*cm);
    
    // hadron calorimeter cell
    G4VSolid* HadCalCellSolid
      = new G4Box("HadCalCellBox",15.*cm,15.*cm,50.*cm);
    G4LogicalVolume* HadCalCellLogical
      = new G4LogicalVolume(HadCalCellSolid,lead,"HadCalCellLogical");
    new G4PVReplica("HadCalCellPhysical",HadCalCellLogical,
                    HadCalColumnLogical,kYAxis,2,30.*cm);
    
    // hadron calorimeter layers
    G4VSolid* HadCalLayerSolid
      = new G4Box("HadCalLayerBox",15.*cm,15.*cm,2.5*cm);
    G4LogicalVolume* HadCalLayerLogical
      = new G4LogicalVolume(HadCalLayerSolid,lead,"HadCalLayerLogical");
    new G4PVReplica("HadCalLayerPhysical",HadCalLayerLogical,
                    HadCalCellLogical,kZAxis,20,5.*cm);
    
    // scintillator plates
    G4VSolid* HadCalScintiSolid
      = new G4Box("HadCalScintiBox",15.*cm,15.*cm,0.5*cm);
    fHadCalScintiLogical
      = new G4LogicalVolume(HadCalScintiSolid,scintillator,
                            "HadCalScintiLogical");
    new G4PVPlacement(0,G4ThreeVector(0.,0.,2.*cm),fHadCalScintiLogical,
                      "HadCalScintiPhysical",HadCalLayerLogical,
                      false,0,checkOverlaps);
    
    // visualization attributes ------------------------------------------------
    
    G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
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
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    
    G4VSensitiveDetector* hodoscope1 
      = new B5HodoscopeSD(SDname="/hodoscope1");
    SDman->AddNewDetector(hodoscope1);
    fHodoscope1Logical->SetSensitiveDetector(hodoscope1);

    G4VSensitiveDetector* hodoscope2 
      = new B5HodoscopeSD(SDname="/hodoscope2");
    SDman->AddNewDetector(hodoscope2);
    fHodoscope2Logical->SetSensitiveDetector(hodoscope2);
    
    G4VSensitiveDetector* chamber1 
      = new B5DriftChamberSD(SDname="/chamber1");
    SDman->AddNewDetector(chamber1);
    fWirePlane1Logical->SetSensitiveDetector(chamber1);

    G4VSensitiveDetector* chamber2 
      = new B5DriftChamberSD(SDname="/chamber2");
    SDman->AddNewDetector(chamber2);
    fWirePlane2Logical->SetSensitiveDetector(chamber2);
    
    G4VSensitiveDetector* emCalorimeter 
      = new B5EmCalorimeterSD(SDname="/EMcalorimeter");
    SDman->AddNewDetector(emCalorimeter);
    fCellLogical->SetSensitiveDetector(emCalorimeter);
    
    G4VSensitiveDetector* hadCalorimeter 
      = new B5HadCalorimeterSD(SDname="/HadCalorimeter");
    SDman->AddNewDetector(hadCalorimeter);
    fHadCalScintiLogical->SetSensitiveDetector(hadCalorimeter);

    // magnetic field ----------------------------------------------------------
    fMagneticField = new B5MagneticField();
    fFieldMgr = new G4FieldManager();
    fFieldMgr->SetDetectorField(fMagneticField);
    fFieldMgr->CreateChordFinder(fMagneticField);
    G4bool forceToAllDaughters = true;
    fMagneticLogical->SetFieldManager(fFieldMgr, forceToAllDaughters);

    // Register the field and its manager for deleting
    G4AutoDelete::Register(fMagneticField);
    G4AutoDelete::Register(fFieldMgr);
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::ConstructMaterials()
{
    G4NistManager* nistManager = G4NistManager::Instance();

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
    // G4Material* air = G4Material::GetMaterial("G4_AIR");
    // G4double density = 1.0e-5*air->GetDensity();
    // nistManager
    //   ->BuildMaterialWithNewDensity("Air_lowDensity", "G4_AIR", density);

    G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::SetArmAngle(G4double val)
{
    if (!fSecondArmPhys)
    {
        G4cerr << "Detector has not yet been constructed." << G4endl;
        return;
    }
    
    fArmAngle = val;
    *fArmRotation = G4RotationMatrix();  // make it unit vector
    fArmRotation->rotateY(fArmAngle);
    G4double x = -5.*m * std::sin(fArmAngle);
    G4double z = 5.*m * std::cos(fArmAngle);
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
    G4GenericMessenger::Command& armAngleCmd
      = fMessenger->DeclareMethodWithUnit("armAngle","deg",
                                  &B5DetectorConstruction::SetArmAngle, 
                                  "Set rotation angle of the second arm.");
    armAngleCmd.SetParameterName("angle", true);
    armAngleCmd.SetRange("angle>=0. && angle<180.");
    armAngleCmd.SetDefaultValue("30.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
