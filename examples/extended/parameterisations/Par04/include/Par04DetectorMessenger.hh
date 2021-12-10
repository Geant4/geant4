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

#ifndef PAR04DETECTORMESSENGER_H
#define PAR04DETECTORMESSENGER_H

#include <G4String.hh>       // for G4String
#include "G4UImessenger.hh"  // for G4UImessenger
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcommand;
class G4UIdirectory;
class Par04DetectorConstruction;

/**
 * @brief Detector messenger.
 *
 * Provides UI commands to setup detector and readout geometry (prior to
 * initialization). Radius, length, and material of the detector, as well as
 * segmentation of the readout geometry can be changed.
 *
 */

class Par04DetectorMessenger : public G4UImessenger
{
 public:
  Par04DetectorMessenger(Par04DetectorConstruction*);
  ~Par04DetectorMessenger();

  /// Invokes appropriate methods based on the typed command
  virtual void SetNewValue(G4UIcommand*, G4String) final;
  /// Retrieves the current settings
  virtual G4String GetCurrentValue(G4UIcommand*) final;

 private:
  /// Detector construction to setup
  Par04DetectorConstruction* fDetector = nullptr;
  /// Command to set the directory common to all messengers in this example
  /// /Par04
  G4UIdirectory* fExampleDir = nullptr;
  /// Command to set the directory for detector settings /Par04/detector
  G4UIdirectory* fDetectorDir = nullptr;
  /// Command printing current settings
  G4UIcmdWithoutParameter* fPrintCmd;
  /// Command to set the detector inner radius
  G4UIcmdWithADoubleAndUnit* fDetectorInnerRadiusCmd = nullptr;
  /// Command to set the detector length
  G4UIcmdWithADoubleAndUnit* fDetectorLengthCmd = nullptr;
  /// Command to set the number of layers
  G4UIcmdWithAnInteger* fNbLayersCmd = nullptr;
  ///  Commanbd to set the absorbers within layers (material, thickness, sensitivity)
  G4UIcommand* fAbsorCmd = nullptr;
  /// Command to set the directory for sensitive detector settings /Par04/mesh
  G4UIdirectory* fMeshDir = nullptr;
  /// Command to set the number of cells in the cylindrical readout mesh (along rho axis)
  G4UIcmdWithAnInteger* fMeshNbRhoCellsCmd = nullptr;
  /// Command to set the number of cells in the cylindrical readout mesh (along phi axis)
  G4UIcmdWithAnInteger* fMeshNbPhiCellsCmd = nullptr;
  /// Command to set the number of cells in the cylindrical readout mesh (along z axis)
  G4UIcmdWithAnInteger* fMeshNbZCellsCmd = nullptr;
  /// Command to the size of cells in the cylindrical readout mesh (along rho axis)
  G4UIcmdWithADoubleAndUnit* fMeshSizeRhoCellsCmd = nullptr;
  /// Command to the size of cells in the cylindrical readout mesh (along z axis)
  G4UIcmdWithADoubleAndUnit* fMeshSizeZCellsCmd = nullptr;
};

#endif
