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

#ifdef USE_INFERENCE_TORCH
#ifndef PAR04TORCHINFERENCE_HH
#define PAR04TORCHINFERENCE_HH
#include <G4String.hh>                         // for G4String
#include <G4Types.hh>                          // for G4int, G4double
#include <memory>                              // for unique_ptr
#include <vector>                              // for vector
#include "Par04InferenceInterface.hh"          // for Par04InferenceInterface
#include <torch/script.h>

/**
 * @brief Inference using the TORCH.
 *
 * Runs the inference with LibTorch using the input vector from Par04InferenceSetup.
 *
 **/

class Par04TorchInference : public Par04InferenceInterface
{
public:
  Par04TorchInference(G4String);
  Par04TorchInference();

  /// Run inference
  /// @param[in] aGenVector Input latent space and conditions
  /// @param[out] aEnergies Model output = generated shower energies
  /// @param[in] aSize Size of the output
  void RunInference(std::vector<float> aGenVector, std::vector<G4double>& aEnergies, int aSize);

private:
  torch::jit::script::Module fModule;

};

#endif /* PAR04TORCHINFERENCE_HH */
#endif
