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

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4CrystalExtension.hh"
#include "G4ExtendedMaterial.hh"
#include "G4LogicalCrystalVolume.hh"

#include "G4ChannelingMaterialData.hh"
#include "G4ChannelingOptrMultiParticleChangeCrossSection.hh"

#include "SensitiveDetector.hh"

#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():
fECfileName("Si220pl"),
fMaterialName("G4_Si"),
fSizes(G4ThreeVector(1.*CLHEP::mm,
                     70.*CLHEP::mm,
                     1.94 * CLHEP::mm)),
fBR(G4ThreeVector(0.,0.,0.)),
fAngles(G4ThreeVector(0.,0.,0.)){
    fMessenger = new DetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::~DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::DefineMaterials(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct(){
    
    //** World **//
    G4Material* worldMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    
    G4double worldSizeXY =  1. * CLHEP::meter;
    G4double worldSizeZ  = 22. * CLHEP::meter;
    
    G4Box* worldSolid = new G4Box("world.solid",
                                  worldSizeXY/2.,
                                  worldSizeXY/2.,
                                  worldSizeZ/2.);
    
    G4LogicalVolume* worldLogic = new G4LogicalVolume(worldSolid,
                                                      worldMaterial,
                                                      "world.logic");
    
    G4PVPlacement* worldPhysical = new G4PVPlacement(0,
                                                     G4ThreeVector(),
                                                     worldLogic,
                                                     "world.physic",
                                                     0,
                                                     false,
                                                     0);
    
    
    //** Detectors instantiation **//
    G4ThreeVector fDetectorSizes(G4ThreeVector(38.0 * CLHEP::mm,
                                               38.0 * CLHEP::mm,
                                               0.64 * CLHEP::mm));
    G4double fDetectorDistance[3] = {
        -9.998 * CLHEP::m,
        -0.320 * CLHEP::m,
        +10.756* CLHEP::m};
                                 
    G4Box* ssdSolid = new G4Box("ssd.solid",
                                fDetectorSizes.x()/2.,
                                fDetectorSizes.y()/2.,
                                fDetectorSizes.z()/2.);
    
    G4Material* detectorMaterial =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    G4LogicalVolume* ssdLogic =
      new G4LogicalVolume(ssdSolid,
                          detectorMaterial,
                          "ssd.logic");
    
    for(size_t i1=0;i1<3;i1++){
        new G4PVPlacement(0,
                          G4ThreeVector(0.,
                                        0.,
                                        fDetectorDistance[i1]),
                          ssdLogic,
                          "ssd.physic",
                          worldLogic,
                          false,
                          i1);
        
    }
    
    //** Crystal solid parameters **//
    G4Box* crystalSolid = new G4Box("crystal.solid",
                                    fSizes.x()/2.,
                                    fSizes.y()/2.,
                                    fSizes.z()/2.);
    
    //** Crystal Definition Start **//
    G4Material* mat =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
    G4ExtendedMaterial* Crystal =
      new G4ExtendedMaterial("crystal.material",mat);
    
    Crystal->RegisterExtension(std::unique_ptr<G4CrystalExtension>(
            new G4CrystalExtension(Crystal)));
    G4CrystalExtension* crystalExtension =
      (G4CrystalExtension*)Crystal->RetrieveExtension("crystal");
    crystalExtension->SetUnitCell(
        new G4CrystalUnitCell(5.43 * CLHEP::angstrom,
                              5.43 * CLHEP::angstrom,
                              5.43 * CLHEP::angstrom,
                              CLHEP::halfpi,
                              CLHEP::halfpi,
                              CLHEP::halfpi,
                              227));
    
    Crystal->RegisterExtension(std::unique_ptr<G4ChannelingMaterialData>(
            new G4ChannelingMaterialData("channeling")));
    G4ChannelingMaterialData* crystalChannelingData =
      (G4ChannelingMaterialData*)Crystal->RetrieveExtension("channeling");
    crystalChannelingData->SetFilename(fECfileName);

    if(fBR!=G4ThreeVector()){
        crystalChannelingData->SetBR(fBR.x());
    }
    
    G4LogicalCrystalVolume* crystalLogic =
      new G4LogicalCrystalVolume(crystalSolid,
                                 Crystal,
                                 "crystal.logic");
    crystalLogic->SetVerbose(1);
    //** Crystal Definition End **//
    
    
    G4RotationMatrix* rot = new G4RotationMatrix;
    if(fAngles.x()!=0.){
        rot->rotateX(fAngles.x());
    }
    if(fAngles.y()!=0.){
        rot->rotateY(fAngles.y());
    }
    if(fAngles.z()!=0.){
        rot->rotateZ(fAngles.z());
    }
    
    new G4PVPlacement(rot,
                      G4ThreeVector(),
                      crystalLogic,
                      "crystal.physic",
                      worldLogic,
                      false,
                      0);
    
#ifndef G4MULTITHREADED
    G4ChannelingOptrMultiParticleChangeCrossSection* testMany =
    new G4ChannelingOptrMultiParticleChangeCrossSection();
    testMany->AttachTo(crystalLogic);
    G4cout << " Attaching biasing operator " << testMany->GetName()
    << " to logical volume " << crystalLogic->GetName()
    << G4endl;
    
    G4VSensitiveDetector* telescope = new SensitiveDetector("/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    for(unsigned int i1=0;i1<3;i1++){
        ssdLogic->SetSensitiveDetector(telescope);
    }
#endif

    return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
void DetectorConstruction::ConstructSDandField(){
    G4LogicalVolume* crystalLogic =
      G4LogicalVolumeStore::GetInstance()->GetVolume("crystal.logic");
    G4ChannelingOptrMultiParticleChangeCrossSection* testMany =
    new G4ChannelingOptrMultiParticleChangeCrossSection();
    testMany->AttachTo(crystalLogic);
    G4cout << " Attaching biasing operator " << testMany->GetName()
    << " to logical volume " << crystalLogic->GetName()
    << G4endl;
    
    G4LogicalVolume* ssdLogic =
      G4LogicalVolumeStore::GetInstance()->GetVolume("ssd.logic");
    G4VSensitiveDetector* telescope = new SensitiveDetector("/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    for(unsigned int i1=0;i1<3;i1++){
        ssdLogic->SetSensitiveDetector(telescope);
    }
}
#else
void DetectorConstruction::ConstructSDandField(){;}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

