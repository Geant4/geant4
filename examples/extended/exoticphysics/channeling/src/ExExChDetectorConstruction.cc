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
//

#include "ExExChDetectorConstruction.hh"
#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"
#include "ExExChSensitiveDetector.hh"

#include "XLatticeManager3.hh"

#include "XLogicalAtomicLattice.hh"
#include "XLogicalAtomicLatticeDiamond.hh"
#include "XLogicalBase.hh"
#include "XUnitCell.hh"

#include "XCrystalPlanarMolierePotential.hh"
#include "XCrystalPlanarMoliereElectricField.hh"
#include "XCrystalPlanarNucleiDensity.hh"
#include "XCrystalPlanarMoliereElectronDensity.hh"
#include "XCrystalPlanarMoliereTempPotential.hh"
#include "XCrystalPlanarMoliereElectricField.hh"
#include "XCrystalPlanarNucleiDensity.hh"
#include "XCrystalPlanarMoliereElectronDensity.hh"
#include "XCrystalCharacteristicArray.hh"
#include "XCrystalIntegratedDensityPlanar.hh"
#include "XCrystalIntegratedDensityHub.hh"
#include "XCrystalIntegratedDensityHub.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExExChDetectorConstruction::ExExChDetectorConstruction():
fWorldLogic(0),fXtalLogic(0){
    
    bXtal = true;
    fXtalAngle = G4ThreeVector(0.,0.,0.);
    fXtalSize = G4ThreeVector(1. * CLHEP::millimeter,
                              70. * CLHEP::millimeter,
                              1.94 * CLHEP::millimeter);
    fXtalCurvatureRadius = G4ThreeVector(38.416 * CLHEP::meter,
                                         0. * CLHEP::meter,
                                         0. * CLHEP::meter);
    
    fXtalCellSize = G4ThreeVector(5.431 * CLHEP::angstrom,
                                  5.431 * CLHEP::angstrom,
                                  5.431 * CLHEP::angstrom);
    
    fXtalCellAngle = G4ThreeVector(90.*CLHEP::deg,
                                   90.*CLHEP::deg,
                                   90.*CLHEP::deg);
    
    fXtalTVA = 0.075 * CLHEP::angstrom;
    fXtalMiller = G4ThreeVector(2,2,0);
    
    SetXtalMaterial("G4_Si");
    
    fMessenger = new ExExChDetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ExExChDetectorConstruction::~ExExChDetectorConstruction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ExExChDetectorConstruction::DefineMaterials(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExExChDetectorConstruction::Construct(){
    //** World **//
    fWorldSizeXY = 20. * CLHEP::centimeter;
    fWorldSizeZ = 2.2 * CLHEP::meter;
    fWorldMaterial = G4NistManager::
        Instance()->FindOrBuildMaterial("G4_Galactic");
    
    fWorldSolid = new G4Box("World",
                            fWorldSizeXY/2.,
                            fWorldSizeXY/2.,
                            fWorldSizeZ/2.);
    
    fWorldLogic = new G4LogicalVolume(fWorldSolid,
                                      fWorldMaterial,
                                      "World");
    
    fWorldPhysical = new G4PVPlacement(0,
                                       G4ThreeVector(),
                                       fWorldLogic,
                                       "World",
                                       0,
                                       false,
                                       0);
    
    //** SSD **//
    G4Material* Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    
    fSSDSize = G4ThreeVector(1.92 * CLHEP::centimeter,
                             1.92 * CLHEP::centimeter,
                             0.06 * CLHEP::centimeter); //
    fSSD0XtalDistance = - 1.0 * CLHEP::meter;
    fSSD1XtalDistance = - 0.01 * CLHEP::meter;
    fSSD2XtalDistance = + 1.0 * CLHEP::meter;
    
    fSSDSolid = new G4Box("SiSD",
                          fSSDSize.x()/2.,
                          fSSDSize.y()/2.,
                          fSSDSize.z()/2.);
    
    fSSDLogic = new G4LogicalVolume(fSSDSolid,Si,"SiSD");
    
    new G4PVPlacement(0,
                      G4ThreeVector(0.,0.,fSSD0XtalDistance),
                      fSSDLogic,"SiSD",
                      fWorldLogic,
                      false,
                      0);
    
    new G4PVPlacement(0,
                      G4ThreeVector(0.,0.,fSSD1XtalDistance),
                      fSSDLogic,
                      "SiSD",
                      fWorldLogic,
                      false,
                      1);
    
    new G4PVPlacement(0,
                      G4ThreeVector(0.,0.,fSSD2XtalDistance),
                      fSSDLogic,
                      "SiSD",
                      fWorldLogic,
                      false,2
                      );
    
#ifndef G4MULTITHREADED
    G4String SDname;
    G4VSensitiveDetector* telescope =
        new ExExChSensitiveDetector(SDname="/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    fSSDLogic->SetSensitiveDetector(telescope);
#endif
    
    if(bXtal) ConstructXtalTarget();
    
    return fWorldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
void ExExChDetectorConstruction::ConstructSDandField(){
    G4String SDname;
    G4VSensitiveDetector* telescope =
        new ExExChSensitiveDetector(SDname="/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    fSSDLogic->SetSensitiveDetector(telescope);
}
#else
void ExExChDetectorConstruction::ConstructSDandField(){
}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::ConstructXtalTarget(){
    if(fXtalCurvatureRadius.x() != 0.){
        double fXtalAngleOut =
            fXtalAngle.y() + fXtalSize.z()/fXtalCurvatureRadius.x();
        fXtalSolid = new G4Tubs("Target",
                                fXtalCurvatureRadius.x() - fXtalSize.x()/2,
                                fXtalCurvatureRadius.x() + fXtalSize.x()/2,
                                fXtalSize.y()/2.,
                                fXtalAngle.y(),
                                fXtalAngleOut);
    }
    else{
        fXtalSolid = new G4Box("Target",
                               fXtalSize.x()/2.,
                               fXtalSize.y()/2.,
                               fXtalSize.z()/2.);
    }
    
    fXtalLogic = new G4LogicalVolume(fXtalSolid,fXtalMaterial,"Target");
    
    G4RotationMatrix* vRotationMatrix =
        new G4RotationMatrix; // Rotates X and Z axes only
    
    if(fXtalCurvatureRadius.x() != 0.){
        vRotationMatrix->rotateX(fXtalAngle.x()-CLHEP::pi*0.5);
        vRotationMatrix->rotateY(fXtalAngle.y());
        vRotationMatrix->rotateZ(fXtalAngle.z());
        fXtalPhysical =
            new G4PVPlacement(vRotationMatrix,
                              G4ThreeVector(-fXtalCurvatureRadius.x(),0.,0.),
                              fXtalLogic,"Target",
                              fWorldLogic,
                              false,
                              0);
    }
    else{
        vRotationMatrix->rotateX(fXtalAngle.x());
        vRotationMatrix->rotateY(fXtalAngle.y());
        vRotationMatrix->rotateZ(fXtalAngle.z());
        fXtalPhysical =
            new G4PVPlacement(vRotationMatrix,
                              G4ThreeVector(0.,0.,0.),
                              fXtalLogic,"Target",
                              fWorldLogic,
                              false,
                              0);
    }
    
    //----------------------------------------
    // Create XLogicalLattice
    //----------------------------------------
    XLogicalLattice* logicalLattice = new XLogicalLattice();
    double vScatteringConstant =
        3.67e-41*CLHEP::second*CLHEP::second*CLHEP::second;
    logicalLattice->SetScatteringConstant(vScatteringConstant);
    
    //----------------------------------------
    // Create XLogicalBase
    //----------------------------------------
    XLogicalAtomicLatticeDiamond *diamond_lattice =
        new XLogicalAtomicLatticeDiamond();
    G4Element* element = G4NistManager::
        Instance()->FindOrBuildElement(G4int(fXtalMaterial->GetZ()));
    XLogicalBase *base = new XLogicalBase(element,diamond_lattice);
    
    //----------------------------------------
    // Create XUnitCell
    //----------------------------------------
    XUnitCell* myCell = new XUnitCell();
    myCell->SetSize(fXtalCellSize);
    myCell->AddBase(base);
    
    //----------------------------------------
    // Create XPhysicalLattice
    //----------------------------------------
    XPhysicalLattice* physicalLattice =
        new XPhysicalLattice(fXtalPhysical, logicalLattice);
    physicalLattice->SetUnitCell(myCell);
    physicalLattice->SetMillerOrientation(G4int(fXtalMiller.x()),
                                          G4int(fXtalMiller.y()),
                                          G4int(fXtalMiller.z()));
    physicalLattice->SetLatticeOrientation(fXtalAngle.x(),
                                           fXtalAngle.y(),
                                           fXtalAngle.z());
    physicalLattice->SetThermalVibrationAmplitude(fXtalTVA);
    physicalLattice->SetCurvatureRadius(fXtalCurvatureRadius);
    
    //----------------------------------------
    // Register XPhysicalLattice
    //----------------------------------------
    if(XLatticeManager3::
        GetXLatticeManager()->GetXPhysicalLattice(fXtalPhysical)
        != physicalLattice){
            XLatticeManager3::
                GetXLatticeManager()->RegisterLattice(physicalLattice);
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.

void ExExChDetectorConstruction::SetXtalMaterial(const G4String& name){
    G4Material* vMaterial = G4Material::GetMaterial(name, false);
    
    if(!vMaterial){
        vMaterial = G4NistManager::Instance()->FindOrBuildMaterial(name);
    }
    
    if (vMaterial && vMaterial != fXtalMaterial) {
        G4cout << "DetectorConstructor::SetXtalMaterial() - New Xtal Material: "
            << vMaterial->GetName() << G4endl;
        fXtalMaterial = vMaterial;
        if(fXtalLogic){
            fXtalLogic->SetMaterial(vMaterial);
            G4RunManager::GetRunManager()->PhysicsHasBeenModified();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String ExExChDetectorConstruction::GetXtalMaterial(){
    if(fXtalMaterial) {
        return fXtalMaterial->GetName();
    }
    return "";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalCurvatureRadius(G4ThreeVector cr){
    if(fXtalCurvatureRadius != cr) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalCurvatureRadius = cr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalSize(G4ThreeVector size) {
    if(fXtalSize != size) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalSize = size;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalAngle(G4ThreeVector angle) {
    if(fXtalAngle != angle) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalAngle = angle;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalCellSize(G4ThreeVector cellsize) {
    if(fXtalCellSize != cellsize) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalCellSize = cellsize;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalMiller(G4ThreeVector miller) {
    if(fXtalMiller != miller) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalMiller = miller;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalCellAngle(
                                            G4ThreeVector cellangle) {
    if(fXtalCellAngle != cellangle) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalCellAngle = cellangle;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExExChDetectorConstruction::SetXtalThermalVibrationAmplitude(
                                            G4double thermvibr) {
    if(fXtalTVA != thermvibr) {
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
        fXtalTVA = thermvibr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..
