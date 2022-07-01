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
#ifndef PAR04MLFASTSIMMODEL_HH
#define PAR04MLFASTSIMMODEL_HH

#include <G4String.hh>                // for G4String
#include <G4ThreeVector.hh>           // for G4ThreeVector
#include <G4Types.hh>                 // for G4bool, G4double
#include <memory>                     // for unique_ptr
#include <vector>                     // for vector
#include "G4VFastSimulationModel.hh"  // for G4VFastSimulationModel
class G4FastSimHitMaker;
class G4FastStep;
class G4FastTrack;
class G4ParticleDefinition;
class G4Region;
class Par04InferenceSetup;

/**
 * @brief Inference for the fast simulation model.
 *
 * Parametrisation of electrons, positrons, and gammas. Energy is deposited and
 * distributed according to ML inference. This class is only a shell that
 * triggers the inference, asks for values, and deposits energies at
 * given positions.
 *
 **/

class Par04MLFastSimModel : public G4VFastSimulationModel
{
 public:
  Par04MLFastSimModel(G4String, G4Region*);
  Par04MLFastSimModel(G4String);
  ~Par04MLFastSimModel();
  /// There are no kinematics constraints. True is returned.
  virtual G4bool ModelTrigger(const G4FastTrack&) final;
  /// Model is applicable to electrons, positrons, and photons.
  virtual G4bool IsApplicable(const G4ParticleDefinition&) final;
  /// Take particle out of the full simulation (kill it at the entrance
  /// depositing all the energy). Calculate energy deposited in the detector
  /// from the NN model inference.
  virtual void DoIt(const G4FastTrack&, G4FastStep&) final;

 private:
  /// Inference model that is NN aware
  Par04InferenceSetup* fInference;
  /// Inference model that is NN aware
  /// Helper class for creation of hits within the sensitive detector
  std::unique_ptr<G4FastSimHitMaker> fHitMaker;
  /// Vector of energy values
  std::vector<G4double> fEnergies;
  /// Vector of positions corresponding to energy values (const for one NN
  /// model)
  std::vector<G4ThreeVector> fPositions;
};
#endif /* PAR04INFERENCEMODEL_HH */
#endif
