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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class


#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4EllipticalTube.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4PVParameterised.hh"

#include "PhysGeoImport.hh"
#include "VoxelParameterisation.hh"
#include "ChemGeoImport.hh"

#include <set>

enum class RunningMode{Phys,Chem};
extern RunningMode gRunMode;

class DetectorConstructionMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

    DetectorConstruction(G4double factor=1, G4int verbose=0, G4bool isVisu=false);

    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;

    G4int GetVerbose(){return fVerbose;}
    void SetVerbose(G4int verbose){fVerbose=verbose;}
    void SetCellDefFilePath(const G4String finput);
    void AddVoxelDefFile(const G4String finput);
    void SetWorldBoxSizes(G4ThreeVector);
    void ParseGeoFileForChemMode(const G4String fn);
    void InsertMoleculeInWorld();
private:

    G4double fFactor{1};
    G4int fVerbose{0};
    G4bool fBVisu{false};

    G4Box* fSolidWorld{nullptr};
    G4LogicalVolume* fLogicWorld=nullptr;
    G4VPhysicalVolume* fPhysWorld=nullptr;

    G4VPhysicalVolume *ConstructFullCellNucleusGeo();
    G4VPhysicalVolume *ConstructVoxelGeo(); // for runing chem separately
    DetectorConstructionMessenger *fDetectorMessenger{nullptr};
    std::unique_ptr<ChemGeoImport> fChemGeoImport{nullptr};

    G4String fCellDefFilePath="";
    std::set<G4String> fVoxelDefFilesList;
    G4double fWorldBoxSizeX = 0., fWorldBoxSizeY =0., fWorldBoxSizeZ =0.;
    G4double fVoxelHalfSizeXYZ{0};
    G4bool fUsingUserDefinedSizesForWorld = false;

    G4Material *fWater{nullptr};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void DetectorConstruction::InsertMoleculeInWorld()
{
    fChemGeoImport->InsertMoleculeInWorld();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
