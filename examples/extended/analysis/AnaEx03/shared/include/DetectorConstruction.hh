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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "CLHEP/Units/SystemOfUnits.h"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

   DetectorConstruction();
   ~DetectorConstruction() override;

  public:
    G4VPhysicalVolume* Construct() override;

    void SetAbsorberMaterial (const G4String&);
    void SetAbsorberThickness(G4double);

    void SetGapMaterial (const G4String&);
    void SetGapThickness(G4double);

    void SetCalorSizeYZ(G4double);
    void SetNbOfLayers (G4int);

    void PrintCalorParameters();

    G4double GetWorldSizeX()  { return fWorldSizeX; }
    G4double GetWorldSizeYZ() { return fWorldSizeYZ; }

    G4double GetCalorThickness() { return fCalorThickness; }
    G4double GetCalorSizeYZ()    { return fCalorSizeYZ; }

    G4int GetNbOfLayers() { return fNbOfLayers; }

    G4Material* GetAbsorberMaterial()  { return fAbsorberMaterial; }
    G4double    GetAbsorberThickness() { return fAbsorberThickness; }

    G4Material* GetGapMaterial()  { return fGapMaterial; }
    G4double    GetGapThickness() { return fGapThickness; }

    const G4VPhysicalVolume* GetphysiWorld() { return fPhysiWorld; }
    const G4VPhysicalVolume* GetAbsorber()   { return fPhysiAbsorber; }
    const G4VPhysicalVolume* GetGap()        { return fPhysiGap; }

  private:
    void DefineMaterials();
    void ComputeCalorParameters();
    G4VPhysicalVolume* ConstructCalorimeter();

    G4Material* fAbsorberMaterial = nullptr;
    G4Material* fGapMaterial = nullptr;
    G4Material* fDefaultMaterial = nullptr;

    G4int      fNbOfLayers = 10;
    G4double   fAbsorberThickness = 10.*CLHEP::mm;
    G4double   fGapThickness = 5.*CLHEP::mm;
    G4double   fCalorSizeYZ = 10.*CLHEP::cm;

    G4double   fCalorThickness = 0.;    // will be computed
    G4double   fLayerThickness= 0.;     // will be computed
    G4double   fWorldSizeYZ = 0.;       // will be computed;
    G4double   fWorldSizeX = 0.;        // will be computed;

    G4Box*             fSolidWorld = nullptr;    //pointer to the solid World
    G4LogicalVolume*   fLogicWorld = nullptr;    //pointer to the logical World
    G4VPhysicalVolume* fPhysiWorld = nullptr;    //pointer to the physical World

    G4Box*             fSolidCalor = nullptr;    //pointer to the solid Calor
    G4LogicalVolume*   fLogicCalor = nullptr;    //pointer to the logical Calor
    G4VPhysicalVolume* fPhysiCalor = nullptr;    //pointer to the physical Calor

    G4Box*             fSolidLayer = nullptr;    //pointer to the solid Layer
    G4LogicalVolume*   fLogicLayer = nullptr;    //pointer to the logical Layer
    G4VPhysicalVolume* fPhysiLayer = nullptr;    //pointer to the physical Layer

    G4Box*             fSolidAbsorber = nullptr; //pointer to the solid Absorber
    G4LogicalVolume*   fLogicAbsorber = nullptr; //pointer to the logical Absorber
    G4VPhysicalVolume* fPhysiAbsorber = nullptr; //pointer to the physical Absorber

    G4Box*             fSolidGap = nullptr;      //pointer to the solid Gap
    G4LogicalVolume*   fLogicGap = nullptr;      //pointer to the logical Gap
    G4VPhysicalVolume* fPhysiGap = nullptr;      //pointer to the physical Gap

    DetectorMessenger* fDetectorMessenger = nullptr;  //pointer to the Messenger
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
    fLayerThickness = fAbsorberThickness + fGapThickness;
    fCalorThickness = fNbOfLayers*fLayerThickness;

    fWorldSizeX = 1.2*fCalorThickness; fWorldSizeYZ = 1.2*fCalorSizeYZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

