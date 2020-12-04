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
#ifndef PAR03DETECTORCONSTRUCTION_H
#define PAR03DETECTORCONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"

class Par03DetectorMessenger;
class G4LogicalVolume;

/**
 * @brief Detector construction.
 *
 * Creates a cylindrical detector, with cylinder axis along Z-axis. It is placed
 * in the world volume so that its bases are located at z=0 and z=Length.
 * Dimensions of the detector (Radius and Length) and material can be set using
 * the UI commands.
 * Readout geometry of the detector is created, and can be set by UI commands.
 * Cells are created along z-axis, azimuthal angle, and radius (cylindrical
 * segmentation).
 * Sensitive detector Par03SensitiveDetector is attached to the
 * cell volume.
 * Region for the detector is created as an envelope of the fast simulation.
 *
 */

class Par03DetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  Par03DetectorConstruction();
  virtual ~Par03DetectorConstruction();

  virtual G4VPhysicalVolume* Construct() final;
  virtual void ConstructSDandField() final;

  // Set radius of the cylindrical detector
  void SetRadius(G4double aRadius);
  // Get radius of the cylindrical detector
  inline G4double GetRadius() const { return fDetectorRadius; };
  // Set length of the cylindrical detector (along z-axis)
  void SetLength(G4double aLength);
  // Get length of the cylindrical detector (along z-axis)
  inline G4double GetLength() const { return fDetectorLength; };
  // Set material of the detector (from NIST materials)
  void SetMaterial(const G4String& aMaterial);
  // Get name of the material of the detector
  inline G4String GetMaterial() const { return fDetectorMaterial->GetName(); };

  // Set number of readout cells along z-axis
  inline void SetNbOfLayers(G4int aNumber) { fNbOfLayers = aNumber; };
  // Get number of readout cells along z-axis
  inline G4int GetNbOfLayers() const { return fNbOfLayers; };
  // Set number of readout cells along radius of cylinder
  inline void SetNbOfRhoCells(G4int aNumber) { fNbOfRhoCells = aNumber; };
  // Get number of readout cells along radius of cylinder
  inline G4int GetNbOfRhoCells() const { return fNbOfRhoCells; };
  // Set number of readout cells in azimuthal angle
  inline void SetNbOfPhiCells(G4int aNumber) { fNbOfPhiCells = aNumber; };
  // Get number of readout cells in azimuthal angle
  inline G4int GetNbOfPhiCells() const { return fNbOfPhiCells; };

  // Print detector information
  void Print() const;

 private:
  /// Messenger that allows to modify geometry
  Par03DetectorMessenger* fDetectorMessenger;
  /// Logical volume of replicated cell
  G4LogicalVolume* fLogicCell = nullptr;
  /// World size (in each X, Y, Z dimension)
  G4double fWorldSize = 10 * m;
  /// Radius of the cylindrical detector
  G4double fDetectorRadius = 10 * cm;
  /// Length of the cylindrical detector (along z axis)
  G4double fDetectorLength = 30 * cm;
  /// Material of the detector
  G4Material* fDetectorMaterial = nullptr;
  /// Number of layers = slices along z axis
  G4int fNbOfLayers = 10;
  /// Number of cells along radius
  G4int fNbOfRhoCells = 10;
  /// Number of cells in azimuthal angle
  G4int fNbOfPhiCells = 10;
};

#endif /* PAR03DETECTORCONSTRUCTION_H */
