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

#ifndef eRositaDetectorConstruction_h
#define eRositaDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

class G4Box;
class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;

class eRositaDetectorConstruction : public G4VUserDetectorConstruction {
public:
    explicit eRositaDetectorConstruction();
    
    ~eRositaDetectorConstruction() override;

    auto Construct() -> G4VPhysicalVolume* override;

    auto GetTracker() -> const G4VPhysicalVolume* {
        return trackerPhysicalVolume;        
    };
    
    // auto GetTrackerFullLength() -> G4double { return trackerFullLength; };
    
    // auto GetTargetFullLength() -> G4double { return targetFullLength; };
    
    // auto GetWorldFullLength() -> G4double { return worldFullLength; };

    void ConstructSDandField() override;

    void SetTargetMaterial(G4String materialName);
    
    void SetTrackerMaterial(G4String materialName);
    
    void SetWorldMaterial(G4String materialName);

private:
    // world
    G4Box* worldSolid; // solid envelope
    G4LogicalVolume* worldLogicalVolume; // logical volume
    G4VPhysicalVolume* worldPhysicalVolume; // physical volume
    G4VisAttributes* worldVisualizationStyle; // visualization style

    // target
    G4Box* targetSolid; // solid
    G4LogicalVolume* targetLogicalVolume; // logical volume
    G4VPhysicalVolume* targetPhysicalVolume; // physical volume
    G4VisAttributes* targetVisualizationStyle; // visualization style

    // tracker
    G4Box* trackerSolid; // solid
    G4LogicalVolume* trackerLogicalVolume; // logical volume
    G4VPhysicalVolume* trackerPhysicalVolume; // physical volume
    G4VisAttributes* trackerVisualizationStyle; // visualization style

    // material
    G4Material* vacuum;
    G4Material* targetMaterial; // target material
    G4Material* trackerMaterial; // tracker material
    G4Material* worldMaterial; // tracker material
    
    // size
    G4double worldHalfLength; // half length of the world volume
    G4double targetHalfLength; // half length of target
    G4double targetHalfDepth; // half depth of target
    G4double trackerHalfLength; // half length of tracker
    G4double trackerHalfDepth; // half depth of tracker

    // spatial position of the target
    G4double targetPositionX; // x coordinate
    G4double targetPositionY; // y coordinate
    G4double targetPositionZ; // z coordinate
    
    // spatial position of the tracker
    G4double trackerPositionX; // x coordinate
    G4double trackerPositionY; // y coordinate
    G4double trackerPositionZ; // z coordinate
};
#endif
