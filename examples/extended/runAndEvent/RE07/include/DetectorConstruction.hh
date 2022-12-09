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
/// \file include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4Cache.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <memory>

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

class G4GlobalMagFieldMessenger;

const G4int kMaxAbsor = 10;  // 0 + 9

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  DetectorConstruction();

 public:
  void SetNbOfAbsor(G4int);
  void SetAbsorMaterial(G4int, const G4String&);
  void SetAbsorThickness(G4int, G4double);

  void SetWorldMaterial(const G4String&);
  void SetCalorSizeYZ(G4double);
  void SetNbOfLayers(G4int);

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

 public:
  void PrintCalorParameters();

  G4double GetWorldSizeX() const { return fWorldSizeX; };
  G4double GetWorldSizeYZ() const { return fWorldSizeYZ; };

  G4double GetCalorThickness() const { return fCalorThickness; };
  G4double GetCalorSizeYZ() const { return fCalorSizeYZ; };

  G4int GetNbOfLayers() const { return fNbOfLayers; };

  G4int GetNbOfAbsor() const { return fNbOfAbsor; };
  G4double GetAbsorThickness(G4int i) const { return fAbsorThickness[i]; };
  const G4Material* GetAbsorMaterial(G4int i) const { return fAbsorMaterial[i]; };

  const G4VPhysicalVolume* GetphysiWorld() const { return fPhysiWorld; };
  const G4Material* GetWorldMaterial() const { return fWorldMaterial; };

 private:
  void ComputeCalorParameters();

  G4int fNbOfAbsor;
  G4Material* fAbsorMaterial[kMaxAbsor];
  G4double fAbsorThickness[kMaxAbsor];

  G4int fNbOfLayers;
  G4double fLayerThickness;

  G4double fCalorSizeYZ;
  G4double fCalorThickness;

  G4Material* fWorldMaterial;
  G4double fWorldSizeYZ;
  G4double fWorldSizeX;

  G4LogicalVolume* fLogicWorld;
  G4VPhysicalVolume* fPhysiWorld;

  G4LogicalVolume* fLogicLayerFront;
  G4LogicalVolume* fLogicLayerBack;

  G4LogicalVolume* fLogicAbsorFront[kMaxAbsor];
  G4LogicalVolume* fLogicAbsorBack[kMaxAbsor];

  std::unique_ptr<DetectorMessenger> fDetectorMessenger;
  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
