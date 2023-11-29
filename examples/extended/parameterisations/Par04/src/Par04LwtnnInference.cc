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
#ifdef USE_INFERENCE_LWTNN
#include "Par04LwtnnInference.hh"
#include <fstream>                     // for ifstream
#include <lwtnn/parse_json.hh>         // for parse_json_graph
#include "Par04InferenceInterface.hh"  // for Par04InferenceInterface

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04LwtnnInference::Par04LwtnnInference(G4String modelPath)
  : Par04InferenceInterface()
{
  // file to read
  std::ifstream input(modelPath);
  // build the graph
  fGraph = std::make_unique<lwt::LightweightGraph>(lwt::parse_json_graph(input));
  input.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04LwtnnInference::RunInference(std::vector<float> aGenVector, std::vector<G4double>& aEnergies,
                                       int aSize)
{
  // generation vector
  fNetworkInputs inputs;
  for(std::size_t i = 0; i < aGenVector.size(); ++i)
  {
    inputs["node_0"]["variable_" + std::to_string(i)] = aGenVector[i];
  }

  // run the inference
  fNetworkOutputs outputs = fGraph->compute(inputs);
  aEnergies.assign(aSize, 0);
  for(int i = 0; i < aSize; i++)
    aEnergies[i] = outputs["out_" + std::to_string(i)];
}

#endif
