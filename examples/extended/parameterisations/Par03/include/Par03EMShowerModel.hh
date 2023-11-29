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
#ifndef PAR03EMSHOWERMODEL_HH
#define PAR03EMSHOWERMODEL_HH

#include "G4VFastSimulationModel.hh"

class Par03EMShowerMessenger;
class G4FastSimHitMaker;

/**
 * @brief Example fast simulation model for EM showers.
 *
 * Parametrisation of electrons, positrons, and gammas. It is triggered if those
 * particles enter the detector so that there is sufficient length for the
 * shower development (max depth, controlled by the UI command).
 *
 * Parametrisation is based on the PDG chapter on the electromagnetic cascades
 * (chapter 33.5). Longitudinal profile of the shower is described with Gamma
 * distribution, with beta parameter on average equal to 0.5 (default value,
 * Fig. 33.21), and alpha parameter calcluated from the incident particle energy
 * and material of the detector (critical energy) following Eq.(33.36).
 *
 * Transverse profile is in this model approximated by the Gaussian
 * distribution, with the mean along the shower axis (incident particle momentum
 * direction) and the standard deviation calculated from the detector material
 * (Moliere radius). This assumes that EM shower is in 90% contained within a
 * cylinder of radius equal to Moliere radius, and that area below Gaussian
 * distribution from `mean-1.645 sigma` to `mean+1.645 sigma` is also equal to
 * 90% of total distribution.
 *
 * Parameters of both distributions (alpha, beta for Gamma, sigma for Gaussian)
 * can be overwritten by UI commands.
 *
 * Parametrisation creates N hits of same energy (N can be set by UI command),
 * using rejection sampling to generate position along shower axis from Gamma
 * distribution, and then sampling from uniform and Gaussian distributions to
 * sample phi and radius, respectively. Created hits are deposited in the
 * detector using its readout geometry, using the helper class G4FastSimHitMaker
 * that locates the volume, and calls appropriate sensitive detector class.
 *
 * PDG Chapter 33:
 * https://pdg.lbl.gov/2019/reviews/rpp2018-rev-passage-particles-matter.pdf
 *
 */

class Par03EMShowerModel : public G4VFastSimulationModel
{
 public:
  Par03EMShowerModel(G4String, G4Region*);
  Par03EMShowerModel(G4String);
  ~Par03EMShowerModel();

  /// There are no kinematics constraints. True is returned.
  virtual G4bool ModelTrigger(const G4FastTrack&) final;
  /// Model is applicable to electrons, positrons, and photons.
  virtual G4bool IsApplicable(const G4ParticleDefinition&) final;
  /// Take particle out of the full simulation (kill it at the entrance
  /// depositing all the energy). Calculate energy deposited in the detector
  /// according to Gamma distribution (along the particle direction) and
  /// Gaussian distribution in the transverse direction. Mean of the Gaussian is
  /// centred on the shower axis. Create energy deposits on a cylindrical mesh.
  /// Parameters of the mesh (size, number of cells) and of the distributions
  /// (alpha, beta for Gamma, sigma for Gaussian) can be set with UI commands.
  virtual void DoIt(const G4FastTrack&, G4FastStep&) final;

  /// Print current settings.
  void Print() const;
  /// Set standard deviation of a Gaussian distribution that describes the
  /// transverse shower profile.
  inline void SetSigma(const G4double aSigma) { fSigma = aSigma; };
  /// Get standard deviation of a Gaussian distribution that describes the
  /// transverse shower profile.
  inline G4double GetSigma() const { return fSigma; };
  /// Set alpha parameter of a Gamma distribution that describes the
  /// longitudinal shower profile.
  inline void SetAlpha(const G4double aAlpha) { fAlpha = aAlpha; };
  /// Get alpha parameter of a Gamma distribution that describes the
  /// longitudinal shower profile.
  inline G4double GetAlpha() const { return fAlpha; };
  /// Set beta parameter of a Gamma distribution that describes the longitudinal
  /// shower profile.
  inline void SetBeta(const G4double aBeta) { fBeta = aBeta; };
  /// Get beta parameter of a Gamma distribution that describes the longitudinal
  /// shower profile.
  inline G4double GetBeta() const { return fBeta; };
  /// Set number of (same energy) hits created in the parametrisation.
  inline void SetNbOfHits(const G4int aNumber) { fNbOfHits = aNumber; };
  /// Get number of (same energy) hits created in the parametrisation.s
  inline G4int GetNbOfHits() const { return fNbOfHits; };
  /// Set maximum depth of shower created in fast simulation. It is expressed in
  /// units of radiaton length.
  inline void SetLongMaxDepth(const G4double aDepth)
  {
    fLongMaxDepth = aDepth;
  };
  /// Get maximum depth of shower created in fast simulation. It is expressed in
  /// units of radiaton length.
  inline G4double GetLongMaxDepth() const { return fLongMaxDepth; };

 private:
  /// Gamma distribution
  inline G4double Gamma(G4double x, G4double alpha, G4double beta)
  {
    return (std::pow(beta, alpha) / std::tgamma(alpha) * std::pow(x, alpha - 1) *
            std::exp(-beta * x));
  }
  /// Gaussian distribution
  inline G4double Gaussian(G4double x, G4double sigma = 1, G4double x0 = 0)
  {
    G4double tmp = (x - x0) / sigma;
    return (1.0 / (std::sqrt(2 * CLHEP::pi) * sigma)) * std::exp(-tmp * tmp / 2);
  }

 private:
  /// Messenger for configuration
  Par03EMShowerMessenger* fMessenger;
  /// Helper class for creation of hits within the sensitive detector
  std::unique_ptr<G4FastSimHitMaker> fHitMaker;
  /// Standard deviation of the Gaussian distribution
  /// Can be changed with UI command `/Par03/fastSim/transverseProfile/sigma
  /// <sigma>`
  /// If sigma is smaller than 0, it will be estimated from the detector
  /// material (Moliere radius).
  G4double fSigma = -1;
  /// Alpha parameter of the Gamma distribution
  /// Can be changed with UI command `/Par03/fastSim/longitudunalProfile/alpha
  /// <alpha>`
  /// If alpha is smaller than 0, it will be estimated from particle energy and
  /// the detector material.
  G4double fAlpha = -1;
  /// Beta parameter of the Gamma distribution
  /// Can be changed with UI command `/Par03/fastSim/longitudinalProfile/beta
  /// <beta>`
  G4double fBeta = 0.5;
  /// Number of (same energy) hits created by the parametrisation. Can be
  /// changed with UI command `/Par03/fastSim/numberOfHits <number>`
  G4int fNbOfHits = 100;
  /// Maximum depth of a shower created in fast simulation.
  /// It is expressed in units of radiation length. Can be changed with UI
  /// command `/Par03/fastSim/longitudinalProfile/maxDepth <depth>`
  G4double fLongMaxDepth = 30;
};
#endif /* PAR03EMSHOWERMODEL_HH */