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
#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4Types.hh"

#include <utility>
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4GenericMessenger;
class HGCalTBMaterials;

/**
 * @brief Detector construction.
 *
 * Creates a detector with a configuration set by UI commands. By default a test
 * module is build (few elements including a silicon sensor).
 * Sensitive detectors are attached to silicon sensors.
 *
 */

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume *Construct();
  virtual void ConstructSDandField();

private:
  /// Define UI commands: configuration of detector (geometry setup), and
  /// maximal step size in silicon
  void DefineCommands();
  /// Set detector setup configuration based on ID
  /// @param[in] aValue ID of detector configuration
  void SelectConfiguration(G4int aValue);
  /// Set maximal size of step within silicon sensor
  /// Can be changed by the UI command /HGCalTestbeam/setup/stepSilicon
  /// @param[in] aValue Max step size
  void SetStepSizeSilicon(G4double aValue);
  /// Helper method placing the logical volumes from map of elements
  void ConstructHGCal();

  /// Pointer to the logical volume of the world
  G4LogicalVolume *fLogicWorld = nullptr;
  /// Pointer to the messenger for UI commands
  G4GenericMessenger *fMessenger = nullptr;
  /// Pointer to the class that defines materials and logical volumes of all
  /// possible elements, to be placed according to the configuration setup
  HGCalTBMaterials *fMaterials = nullptr;
  /// Map of elements to be placed along z axis: name and distance to the
  /// previous element is specified
  std::vector<std::pair<G4String, G4double>> fElementsMap;
  /// Viewpoint for the visualisation
  G4double fVisViewpoint = 0;
  /// Configuration of the detector setup
  /// Can be changed by the UI command /HGCalTestbeam/setup/configuration
  G4int fConfiguration = 0;
};

#endif /* DETECTORCONSTRUCTION_HH */
