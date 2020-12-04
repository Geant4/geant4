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

#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4MoleculeGun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction() 
    : G4VUserDetectorConstruction()
{
    fpDNAParser.reset(new DNAParser());
    fpDNAParser->ParseFile("VoxelStraight.fab2g4dna");
    fpGun = fpDNAParser->ReleaseMoleculeGun();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()
{
//Physics World
    G4NistManager * pMan = G4NistManager::Instance();
    G4Material *pWater = pMan->FindOrBuildMaterial("G4_WATER");
    G4Box* pSolidWorld = new G4Box("solidWorld",200*nm, 200*nm, 200*nm);
    G4LogicalVolume* pLogicWorld = new G4LogicalVolume(pSolidWorld,
                                                       pWater,
                                                       "logicWorld");
    G4VPhysicalVolume* pPhysWorld = new G4PVPlacement(0,
                                                      G4ThreeVector(),
                                                      "physWorld", 
                                                      pLogicWorld,
                                                      0, 
                                                      false, 
                                                      0);

    G4LogicalVolume* pLogicStraightVoxel = fpDNAParser->CreateLogicVolume();
    fGeometryMap = fpDNAParser->ReleaseGeoData();
        
    G4VPhysicalVolume* pPhyVoxel = new G4PVPlacement(0, 
                                                    G4ThreeVector(),
                                                    pLogicStraightVoxel, 
                                                    "VoxelStraight",
                                                    pLogicWorld,
                                                    false,
                                                    0);

    fGeometryMap.insert(std::make_pair(pPhysWorld,
                                       DNAVolumeType::physWorld));
    fGeometryMap.insert(std::make_pair(pPhyVoxel,
                                       DNAVolumeType::VoxelStraight));

    return pPhysWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::UpdateGeometry()
{
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
