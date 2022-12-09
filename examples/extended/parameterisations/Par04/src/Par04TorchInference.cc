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
#include "Par04TorchInference.hh"
#include <algorithm>                           // for copy, max
#include <cassert>                             // for assert
#include <cstddef>                             // for size_t
#include <cstdint>                             // for int64_t
#include <utility>                             // for move
#include "Par04InferenceInterface.hh"          // for Par04InferenceInterface
#include <torch/torch.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04TorchInference::Par04TorchInference(G4String modelPath)
  : Par04InferenceInterface()
{
  fModule = torch::jit::load( modelPath );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04TorchInference::RunInference(std::vector<float> aGenVector, std::vector<G4double>& aEnergies,
                                          int aSize)
{
  // latentSize : size of the latent space
  // 4 is the size of the condition vector
  int latentSize = aGenVector.size() - 4;
  // split into latent and condition vectors
  std::vector<float> latent;
  for ( int i=0;i<latentSize;i++) {
    latent.push_back(aGenVector[i]);
  }
  std::vector<float> energy;
  energy.push_back(aGenVector[latentSize+1]);
  std::vector<float> angle;
  energy.push_back(aGenVector[latentSize+2]);
  std::vector<float> geo;
  for ( int i=latentSize+2;i<latentSize+4;i++) {
    geo.push_back(aGenVector[i]);
  }

  // convert vectors to tensors
  torch::Tensor latentVector = torch::tensor(latent);
  torch::Tensor eTensor = torch::tensor(energy);
  torch::Tensor angleTensor = torch::tensor(angle);
  torch::Tensor geoTensor = torch::tensor(geo);

  std::vector<torch::jit::IValue> genInput;

  genInput.push_back( latentVector );
  genInput.push_back( eTensor );
  genInput.push_back( angleTensor );
  genInput.push_back( geoTensor );

  at::Tensor  outTensor = fModule.forward( genInput).toTensor().contiguous();

  std::vector<G4double> output( outTensor.data_ptr<float>(), outTensor.data_ptr<float>() + outTensor.numel() );

  aEnergies.assign(aSize, 0);
  for(int i = 0; i < aSize; i++) {
    aEnergies[i] = output[i];
  }
}

#endif
