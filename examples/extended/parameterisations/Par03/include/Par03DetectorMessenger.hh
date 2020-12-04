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

#ifndef PAR03DETECTORMESSENGER_H
#define PAR03DETECTORMESSENGER_H

#include "G4UImessenger.hh"

class Par03DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

/**
 * @brief Detector messenger.
 *
 * Provides UI commands to setup detector and readout geometry (prior to
 * initialization). Radius, length, and material of the detector, as well as
 * segmentation of the readout geometry can be changed.
 *
 */

class Par03DetectorMessenger : public G4UImessenger
{
 public:
  Par03DetectorMessenger(Par03DetectorConstruction*);
  ~Par03DetectorMessenger();

  /// Invokes appropriate methods based on the typed command
  virtual void SetNewValue(G4UIcommand*, G4String) final;
  /// Retrieves the current settings
  virtual G4String GetCurrentValue(G4UIcommand*) final;

 private:
  /// Detector construction to setup
  Par03DetectorConstruction* fDetector = nullptr;
  /// Command to set the directory common to all messengers in this example
  /// /Par03
  G4UIdirectory* fExampleDir = nullptr;
  /// Command to set the directory for detector settings /Par03/detector
  G4UIdirectory* fDetectorDir = nullptr;
  /// Command printing current settings
  G4UIcmdWithoutParameter* fPrintCmd;
  /// Command to set the detector radius
  G4UIcmdWithADoubleAndUnit* fDetectorRadiusCmd = nullptr;
  /// Command to set the detector length
  G4UIcmdWithADoubleAndUnit* fDetectorLengthCmd = nullptr;
  /// Command to set the detector material
  G4UIcmdWithAString* fDetectorMaterialCmd = nullptr;
  /// Command to set the number of layers
  G4UIcmdWithAnInteger* fNbLayersCmd = nullptr;
  /// Command to set the number of radial cells
  G4UIcmdWithAnInteger* fNbRhoCellsCmd = nullptr;
  /// Command to set the number of cells in azimuthal angle
  G4UIcmdWithAnInteger* fNbPhiCellsCmd = nullptr;
};

#endif
