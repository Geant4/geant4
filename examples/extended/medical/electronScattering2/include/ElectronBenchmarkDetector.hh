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
/// \file medical/electronScattering2/include/ElectronBenchmarkDetector.hh
/// \brief Definition of the ElectronBenchmarkDetector class

#ifndef ElectronBenchmarkDetector_h
#define ElectronBenchmarkDetector_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Cache.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Tubs;
class G4MultiFunctionalDetector;
class G4Material;
class G4VisAttributes;
class ElectronBenchmarkDetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ElectronBenchmarkDetector : public G4VUserDetectorConstruction
{
public:
    
    ElectronBenchmarkDetector();
    virtual ~ElectronBenchmarkDetector();
    
    virtual G4VPhysicalVolume* Construct();
    
    // Sensitive Detector
    void ConstructSDandField();
    
    // Material Definition
    void DefineMaterials();
    
    // Geometry Definition
    G4VPhysicalVolume* CreateWorld();
    void CreateExitWindow(G4LogicalVolume* logicWorld);
    void CreatePrimaryFoil(G4LogicalVolume* logicWorld);
    void CreateMonitor(G4LogicalVolume* logicWorld);
    void CreateHeliumBag(G4LogicalVolume* logicWorld);
    void CreateScorer(G4LogicalVolume* logicWorld);
    G4VPhysicalVolume* CreateGeometry();
    
    // Command Interface
    void SetPrimFoilMaterial(G4String matname);
    void SetPrimFoilThickness(G4double thicknessPrimFoil);
    
private:
    // Exit Window
    G4double fPosWindow0;
    G4double fPosWindow1;
    
    // Primary Foil
    G4Material* fMaterialPrimFoil;
    G4double fHalfThicknessPrimFoil;
    G4double fPosPrimFoil;
    G4LogicalVolume* fLogPrimFoil;
    G4Tubs* fSolidPrimFoil; 

    // Monitor Chambers
    G4double fPosMon0;
    G4double fPosMon1;
    
    // Helium Bag
    G4double fPosBag0;
    G4double fPosBag1;
    G4double fPosHelium0;
    G4double fPosHelium1;
    G4double fThicknessRing;
    
    // Scoring Plane
    G4double fPosScorer;
    G4double fThicknessScorer;
    G4double fWidthScorerRing;
    G4LogicalVolume* fScorerRingLog;
    
    // Radii
    G4double fRadOverall;
    G4double fRadRingInner;
    
    // Extra space remaining in world volume around apparatus
    G4double fPosDelta;
    G4double fRadDelta;
    
    // World volume
    G4LogicalVolume* fLogWorld;
    
    // SensitiveDetector
    // G4Cache mechanism is necessary for multi-threaded operation
    // as it allows us to store separate detector pointer per thread
    const G4Cache<G4MultiFunctionalDetector*> fSensitiveDetectorCache;
    
    // Messenger
    ElectronBenchmarkDetectorMessenger* fMessenger;
    
    // Visualization Attributes
    G4VisAttributes* fWorldVisAtt;
    G4VisAttributes* fWindowVisAtt;
    G4VisAttributes* fPrimFoilVisAtt;
    G4VisAttributes* fMonVisAtt;
    G4VisAttributes* fBagVisAtt;
    G4VisAttributes* fHeliumVisAtt;
    G4VisAttributes* fRingVisAtt;
    G4VisAttributes* fScorerVisAtt;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
