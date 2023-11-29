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
#ifdef USE_INFERENCE
#ifndef PAR04INFERENCEMESSENGER_H
#define PAR04INFERENCEMESSENGER_H

#include <G4String.hh>           // for G4String
#include "G4UImessenger.hh"      // for G4UImessenger
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcommand;
class G4UIdirectory;
class Par04InferenceSetup;

/**
 * @brief Inference messenger.
 *
 * Provides UI commands to setup the inference: name of the inference library,
 * path to the ML model, size of the latent space, the size of the condition vector and
 * flags for optimization and profiling.
 * It allows to specify the mesh size (in cylindrical coordinates) that was used in training dataset
 * (full simulation).
 *
 */

class Par04InferenceMessenger : public G4UImessenger
{
 public:
  Par04InferenceMessenger(Par04InferenceSetup*);
  ~Par04InferenceMessenger();
  /// Invokes appropriate methods based on the typed command
  virtual void SetNewValue(G4UIcommand*, G4String) final;
  /// Retrieves the current settings
  virtual G4String GetCurrentValue(G4UIcommand*) final;

 private:
  /// Inference to setup
  Par04InferenceSetup* fInference = nullptr;
  /// Command to set the directory common to all inference messengers in this example
  /// /Par04
  G4UIdirectory* fExampleDir = nullptr;
  /// Command to set the directory for inference settings /Par04/inference
  G4UIdirectory* fInferenceDir = nullptr;
  /// Command to set the inference library
  G4UIcmdWithAString* fInferenceLibraryCmd = nullptr;
  /// Command to set fModelPathNameCmd
  G4UIcmdWithAString* fModelPathNameCmd = nullptr;
  /// Command to set the fSizeLatentVectorCmd
  G4UIcmdWithAnInteger* fSizeLatentVectorCmd = nullptr;
  /// Command to set the fSizeConditionVectorCmd
  G4UIcmdWithAnInteger* fSizeConditionVectorCmd = nullptr;
  /// Command to set the fProfileFlagCmd
  G4UIcmdWithAnInteger* fProfileFlagCmd = nullptr;
  /// Command to set the fOptimizationFlagCmd
  G4UIcmdWithAnInteger* fOptimizationFlagCmd = nullptr;
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
#endif
