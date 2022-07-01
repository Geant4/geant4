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
#ifdef USE_INFERENCE_ONNX
#include "Par04OnnxInference.hh"
#include <core/session/onnxruntime_cxx_api.h>  // for Value, Session, Env
#include <algorithm>                           // for copy, max
#include <cassert>                             // for assert
#include <cstddef>                             // for size_t
#include <cstdint>                             // for int64_t
#include <utility>                             // for move
#include "Par04InferenceInterface.hh"          // for Par04InferenceInterface


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04OnnxInference::Par04OnnxInference(G4String modelPath, G4int profileFlag, G4int optimizeFlag,
                                       G4int intraOpNumThreads)
  : Par04InferenceInterface()
{
  // initialization of the enviroment and inference session
  auto envLocal = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ENV");
  fEnv          = std::move(envLocal);
  fSessionOptions.SetIntraOpNumThreads(intraOpNumThreads);
  // graph optimizations of the model
  // if the flag is not set to true none of the optimizations will be applied
  // if it is set to true all the optimizations will be applied
  if(optimizeFlag)
  {
    fSessionOptions.SetOptimizedModelFilePath("opt-graph");
    fSessionOptions.SetGraphOptimizationLevel(ORT_ENABLE_ALL);
    // ORT_ENABLE_BASIC #### ORT_ENABLE_EXTENDED
  }
  else
    fSessionOptions.SetGraphOptimizationLevel(ORT_DISABLE_ALL);
  // save json file for model execution profiling
  if(profileFlag)
    fSessionOptions.EnableProfiling("opt.json");

  auto sessionLocal = std::make_unique<Ort::Session>(*fEnv, modelPath, fSessionOptions);
  fSession          = std::move(sessionLocal);
  fInfo = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemTypeDefault);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04OnnxInference::RunInference(std::vector<float> aGenVector, std::vector<G4double>& aEnergies,
                                      int aSize)
{
  // input nodes
  Ort::AllocatorWithDefaultOptions allocator;
  std::vector<int64_t> input_node_dims;
  size_t num_input_nodes = fSession->GetInputCount();
  std::vector<const char*> input_node_names(num_input_nodes);
  for(std::size_t i = 0; i < num_input_nodes; i++)
  {
    char* input_name               = fSession->GetInputName(i, allocator);
    fInames                        = { input_name };
    input_node_names[i]            = input_name;
    Ort::TypeInfo type_info        = fSession->GetInputTypeInfo(i);
    auto tensor_info               = type_info.GetTensorTypeAndShapeInfo();
    input_node_dims                = tensor_info.GetShape();
    for(std::size_t j = 0; j < input_node_dims.size(); j++)
    {
      if(input_node_dims[j] < 0)
        input_node_dims[j] = 1;
    }
  }
  // output nodes
  std::vector<int64_t> output_node_dims;
  size_t num_output_nodes = fSession->GetOutputCount();
  std::vector<const char*> output_node_names(num_output_nodes);
  for(std::size_t i = 0; i < num_output_nodes; i++)
  {
    char* output_name              = fSession->GetOutputName(i, allocator);
    output_node_names[i]           = output_name;
    Ort::TypeInfo type_info        = fSession->GetOutputTypeInfo(i);
    auto tensor_info               = type_info.GetTensorTypeAndShapeInfo();
    output_node_dims               = tensor_info.GetShape();
    for(std::size_t j = 0; j < output_node_dims.size(); j++)
    {
      if(output_node_dims[j] < 0)
        output_node_dims[j] = 1;
    }
  }

  // create input tensor object from data values
  std::vector<int64_t> dims = { 1, (unsigned) (aGenVector.size()) };
  Ort::Value Input_noise_tensor =
    Ort::Value::CreateTensor<float>(fInfo, aGenVector.data(), aGenVector.size(), dims.data(), dims.size());
  assert(Input_noise_tensor.IsTensor());
  std::vector<Ort::Value> ort_inputs;
  ort_inputs.push_back(std::move(Input_noise_tensor));
  // run the inference session
  std::vector<Ort::Value> ort_outputs =
    fSession->Run(Ort::RunOptions{ nullptr }, fInames.data(), ort_inputs.data(), ort_inputs.size(),
                  output_node_names.data(), output_node_names.size());
  // get pointer to output tensor float values
  float* floatarr = ort_outputs.front().GetTensorMutableData<float>();
  aEnergies.assign(aSize, 0);
  for(int i = 0; i < aSize; ++i)
    aEnergies[i] = floatarr[i];
}

#endif
