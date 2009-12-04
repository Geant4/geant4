//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#ifndef ElectronBenchmarkDetector_h
#define ElectronBenchmarkDetector_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VisAttributes;
class ElectronBenchmarkDetectorMessenger;

class ElectronBenchmarkDetector : public G4VUserDetectorConstruction
{
  public:
  
    ElectronBenchmarkDetector();
   ~ElectronBenchmarkDetector();

    virtual G4VPhysicalVolume* Construct();

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
    void UpdateGeometry();

    // Scorer Activation
    void ActivateScorer();

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
    
    // Messenger
    ElectronBenchmarkDetectorMessenger* fMessenger;
    
    // Visualization Attributes
    G4VisAttributes* worldVisAtt;
    G4VisAttributes* windowVisAtt;
    G4VisAttributes* primFoilVisAtt;
    G4VisAttributes* monVisAtt;
    G4VisAttributes* bagVisAtt;
    G4VisAttributes* heliumVisAtt;
    G4VisAttributes* ringVisAtt;
    G4VisAttributes* scorerVisAtt;
};

#endif
