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
// --------------------------------------------------------------
//

#ifndef ExExChDetectorConstruction_h
#define ExExChDetectorConstruction_h 1
#endif

#include "G4VUserDetectorConstruction.hh"
#include "ExExChDetectorConstructionMessenger.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"

#include "globals.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExExChDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    
    ExExChDetectorConstruction();
    ~ExExChDetectorConstruction();
    
    void DefineMaterials();
    G4VPhysicalVolume* Construct();

private:
    ExExChDetectorConstructionMessenger* fMessenger;
    
    void ConstructSDandField();

private:
    G4ThreeVector fWorldSize;
    G4Box* fWorldSolid;
    G4LogicalVolume* fWorldLogic;
    G4VPhysicalVolume* fWorldPhysical;
    G4Material* fWorldMaterial;

    //** SSD **//
private:
    void ConstructSiliconStripDetectors();
    G4LogicalVolume* ConstructSiSD(G4int);
    G4bool bSiSD;
    G4ThreeVector fSSDSize;
    G4double fSSDXtalDistance[3];
    G4LogicalVolume* fSSDLogic[3];
    G4double fSSDBoxThickness;
    G4ThreeVector fSSDBoxSize;    
    //** Beam pipe **//

    void SetBeamPipeRadius(G4double aDouble) {fBeamPipeRadius = aDouble;};
    G4double GetBeamPipeRadius() {return fBeamPipeRadius;};

    void SetBeamPipeThickness(G4double aDouble) {fBeamPipeThickness = aDouble;};
    G4double GetBeamPipeThickness() {return fBeamPipeThickness;};

private:
    G4LogicalVolume* ConstructBeamPipe(G4double);
    G4bool bBeamPipe;
    G4double fBeamPipeRadius;
    G4double fBeamPipeThickness;
    //** Xtal **//
public:
    void AddXtalTarget() {
        bXtal = true;
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
    };
    void SetXtalMaterial(const G4String& name);
    G4String GetXtalMaterial();
    void SetXtalCurvatureRadius(G4ThreeVector);
    G4ThreeVector GetXtalCurvatureRadius() {return fXtalCurvatureRadius;};
    void SetXtalSize(G4ThreeVector);
    G4ThreeVector GetXtalSize() {return fXtalSize;};
    void SetXtalAngle(G4ThreeVector); 
    G4ThreeVector GetXtalAngle() {return fXtalAngle;};
    void SetXtalCellSize(G4ThreeVector);
    G4ThreeVector GetXtalCellSize() {return fXtalCellSize;};
    void SetXtalCellAngle(G4ThreeVector);
    G4ThreeVector GetXtalCellAngle() {return fXtalCellAngle;};
    void SetXtalThermalVibrationAmplitude(G4double);
    G4double GetXtalThermalVibrationAmplitude() {return fXtalTVA;};
    void SetXtalMiller(G4ThreeVector);
    G4ThreeVector GetXtalMiller() {return fXtalMiller;};
    
private:
    void ConstructXtalTarget();
    G4bool bXtal;
    
    G4ThreeVector fXtalCurvatureRadius;
    
    G4Material* fXtalMaterial;
    G4ThreeVector fXtalAngle;
    G4ThreeVector fXtalSize;
    G4ThreeVector fXtalCellSize;
    G4ThreeVector fXtalCellAngle;
    G4ThreeVector fXtalMiller;
    G4double fXtalTVA;
    
    G4VSolid* fXtalSolid;
    G4LogicalVolume* fXtalLogic;
    G4VPhysicalVolume* fXtalPhysical;
};


