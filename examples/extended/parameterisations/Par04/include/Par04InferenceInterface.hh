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
#ifndef PAR04INFERENCEINTERFACE_HH
#define PAR04INFERENCEINTERFACE_HH

#include "globals.hh"
#include <vector>

/**
 * @brief Inference interface
 *
 * Base class interface that allows to read in ML models and run the inference.
 *
 */

class Par04InferenceInterface
{
 public:
  virtual ~Par04InferenceInterface(){};

  /// Run inference
  /// @param[in] aGenVector Input latent space and conditions
  /// @param[out] aEnergies Model output = generated shower energies
  /// @param[in] aSize Size of the output
  virtual void RunInference(std::vector<float> aGenVector, std::vector<G4double>& aEnergies,
                            int aSize) = 0;
};

#endif /* PAR04INFERENCEINTERFACE_HH */
#endif
