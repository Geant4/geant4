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
#ifndef PAR04DETECTORCONSTRUCTION_H
#define PAR04DETECTORCONSTRUCTION_H

#include <CLHEP/Units/SystemOfUnits.h>     // for cm, mm, pi, rad
#include <G4String.hh>                     // for G4String
#include <G4Types.hh>                      // for G4double, G4bool, G4int
#include <array>                           // for array
#include <cstddef>                         // for size_t
#include <vector>                          // for vector
#include "G4Material.hh"                   // for G4Material
#include "G4SystemOfUnits.hh"              // for cm, mm, rad
#include "G4ThreeVector.hh"                // for G4ThreeVector
#include "G4VUserDetectorConstruction.hh"  // for G4VUserDetectorConstruction
class G4LogicalVolume;
class G4VPhysicalVolume;
class Par04DetectorMessenger;

/**
 * @brief Detector construction.
 *
 * Creates a cylindrical detector, with cylinder axis along Z-axis. It is placed
 * in the  centre of the world volume.
 * Dimensions of the detector (inner radius and its length) as well as composition (number of
 * radial layers, number of absorbers, its thicknesses and materials) can be set using the UI
 * commands. There may be up to two absorbers used to build the layers.
 *
 * TODO Extend to allow use of more than two different absorbers.
 *
 * Each absorber may have differnt thickness (along radial direction), material, and may be
 * either sensitive to particle passage (will register deposited energy) or not (passive).
 * Readout geometry of the detector is dynamic, and its size can be set by UI commands.
 * Cells are created along z-axis, azimuthal angle, and radius (cylindrical
 * segmentation). The z axis is parallel to the direction of the particle entering the
 * detector volume. The mesh also starts at the entrance to the detector volume.
 * In order to define this enrance position and direction, a fast simulation model is used.
 * For the detector volume, if particle has entered, sets the particle direction and position
 * in the event information. Those vectors define the readout mesh for the event.
 *
 * TODO In order to speed up the simulation, fast simulation that checks for the entrance
 * properites should be defined in a very thin region instead of a whole region of the detector.
 *
 * Sensitive detector Par04SensitiveDetector is attached to the detector volume
 * Region for the detector is created as an envelope of the fast simulation.
 *
 */

class Par04DetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  Par04DetectorConstruction();
  virtual ~Par04DetectorConstruction();

  virtual G4VPhysicalVolume* Construct() final;
  virtual void ConstructSDandField() final;

  ///  Set inner radius of the cylindrical detector
  void SetInnerRadius(G4double aInnerRadius);
  /// Get inner radius of the cylindrical detector
  inline G4double GetInnerRadius() const { return fDetectorInnerRadius; };
  ///  Set length radius of the cylindrical detector
  void SetLength(G4double aLength);
  /// Get length of the cylindrical detector (along z-axis)
  inline G4double GetLength() const { return fDetectorLength; };
  ///  Set number of layers
  inline void SetNbOfLayers(G4int aNumber) { fNbOfLayers = aNumber; };
  ///  Get number of layers
  inline G4int GetNbOfLayers() const { return fNbOfLayers; };

  ///  Set material of the layer (from NIST materials)
  void SetAbsorberMaterial(const std::size_t aLayer, const G4String& aMaterial);
  ///  Get name of the material of the layer
  inline G4String GetAbsorberMaterial(const std::size_t aLayer) const
  {
    return fAbsorberMaterial[aLayer]->GetName();
  };
  ///  Set thickness of the layer
  void SetAbsorberThickness(const std::size_t aLayer, const G4double aThickness);
  ///  Get thickness of the layer
  inline G4double GetAbsorberThickness(const std::size_t aLayer) const
  {
    return fAbsorberThickness[aLayer];
  };
  ///  Set sensitivity of the layer
  void SetAbsorberSensitivity(const std::size_t aLayer, const G4bool aSensitivity);
  ///  Get sensitivity of the layer
  inline G4bool GetAbsorberSensitivity(const std::size_t aLayer) const
  {
    return fAbsorberSensitivity[aLayer];
  };

  /// Set number of Mesh cells in cylindrical coordinates (r, phi, z)
  inline void SetMeshNbOfCells(G4ThreeVector aNb) { fMeshNbOfCells = aNb; };
  /// Set number of Mesh cells in cylindrical coordinates along one of the axis
  /// @param[in] aIndex index of cylindrical axis (0,1,2) = (r, phi, z)
  inline void SetMeshNbOfCells(std::size_t aIndex, G4double aNb) { fMeshNbOfCells[aIndex] = aNb; };
  /// Get number of Mesh cells in cylindrical coordinates (r, phi, z)
  inline G4ThreeVector GetMeshNbOfCells() const { return fMeshNbOfCells; };
  /// Set size of Mesh cells in cylindrical coordinates (r, phi, z)
  inline void SetMeshSizeOfCells(G4ThreeVector aNb) { fMeshSizeOfCells = aNb; };
  /// Set size of Mesh cells in cylindrical coordinates along one of the axis
  /// @param[in] aIndex index of cylindrical axis (0,1,2) = (r, phi, z)
  inline void SetMeshSizeOfCells(std::size_t aIndex, G4double aNb) { fMeshSizeOfCells[aIndex] = aNb; };
  /// Get size of Mesh cells in cylindrical coordinates (r, phi, z)
  inline G4ThreeVector GetMeshSizeOfCells() const { return fMeshSizeOfCells; };

  /// Print detector information
  void Print() const;

 private:
  ///  Messenger that allows to modify geometry
  Par04DetectorMessenger* fDetectorMessenger = nullptr;
  ///  Inner radius of the cylindrical detector
  G4double fDetectorInnerRadius = 80 * cm;
  ///  Length of the cylindrical detector (along z axis)
  G4double fDetectorLength = 24 * cm;
  ///  Logical volume(s) of the sensitive absorbers
  std::vector<G4LogicalVolume*> fLayerLogical;
  ///  Material(s) of the layers
  std::array<G4Material*, 2> fAbsorberMaterial = { nullptr, nullptr };
  ///  Thickness(es) of the layers
  std::array<G4double, 2> fAbsorberThickness = { 1 * cm, 0 };
  ///  Sensitivity of the layers
  std::array<G4bool, 2> fAbsorberSensitivity = { true, 0 };
  ///  Number of layers = slices along z axis
  G4int fNbOfLayers = 24;
  ///  Mesh number of cells (Nr, Nphi, Nz)
  G4ThreeVector fMeshNbOfCells = { 40, 50, 48 };
  ///  Mesh size of cells (dr, dphi, dz).
  G4ThreeVector fMeshSizeOfCells = { 5 * mm, 2 * CLHEP::pi / 50 * CLHEP::rad, 5 * mm };
};

#endif /* PAR04DETECTORCONSTRUCTION_H */
