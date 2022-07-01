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
/// \brief Detector Construction for a DNA molecule in a cell

#ifndef MOLECULAR_DETECTOR_CONSTRUCTION_HH
#define MOLECULAR_DETECTOR_CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VDNAMolecularGeometry.hh"

#include <map>

class DetectorMessenger;

class DNAGeometry;

class G4LogicalVolume;

class G4VPhysicalVolume;

class G4Material;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  DetectorConstruction();

  ~DetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;

  void ConstructSDandField() override{};

  void SetWorldSideLength(const G4double& length) { fWorldSideLength = length; }

  void SetCellRadius(const G4ThreeVector& length) { fCellRadius = length; }

  [[maybe_unused]] auto GetWorld() const { return fWorld; };

  auto GetDNAGeometry() const { return fpDNAGeometry; };

 protected:
  G4VPhysicalVolume* ConstructDetector();

  void BuildMaterials();

 private:
  // G4bool fCheckOverlaps;
  G4VPhysicalVolume* fWorld = nullptr;
  DNAGeometry* fpDNAGeometry;

  // Misc messenger values
  DetectorMessenger* fpDetectorMessenger;
  G4double fWorldSideLength = 6 * um;
  G4ThreeVector fCellRadius;

  // materials
  G4Material* fpWorld = nullptr;  // dousatsu
  G4Material* fpWater = nullptr;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_DETECTOR_CONSTRUCTION_HH
