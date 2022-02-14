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
// -------------------------------------------------------------------
//
// Geant4 Class header file
//
// File name:     G4PolarizedIonisationModel
//
// Author:        A.Schaelicke on base of Vladimir Ivanchenko code
//
// Class Description:
//   Physics implementation of polarized Bhabha/Moller scattering
//
// -------------------------------------------------------------------

#ifndef G4PolarizedIonisationModel_h
#define G4PolarizedIonisationModel_h 1

#include "globals.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4StokesVector.hh"
#include "G4ThreeVector.hh"

class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4ParticleDefinition;
class G4VPolarizedXS;

class G4PolarizedIonisationModel : public G4MollerBhabhaModel
{
 public:
  explicit G4PolarizedIonisationModel(
    const G4ParticleDefinition* p = nullptr,
    const G4String& nam           = "PolarizedMollerBhabha");

  virtual ~G4PolarizedIonisationModel() override;

  virtual G4double ComputeCrossSectionPerElectron(const G4ParticleDefinition*,
                                                  G4double kinEnergy,
                                                  G4double cut,
                                                  G4double emax) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*, G4double tmin,
                                 G4double maxEnergy) override;

  G4PolarizedIonisationModel(G4PolarizedIonisationModel&) = delete;
  G4PolarizedIonisationModel& operator=(
    const G4PolarizedIonisationModel& right) = delete;

  void SetTargetPolarization(const G4ThreeVector& pTarget)
  {
    fTargetPolarization = G4StokesVector(pTarget);
  }
  void SetBeamPolarization(const G4ThreeVector& pBeam)
  {
    fBeamPolarization = G4StokesVector(pBeam);
  }
  const G4StokesVector& GetTargetPolarization() { return fTargetPolarization; }
  const G4StokesVector& GetBeamPolarization() { return fBeamPolarization; }
  const G4StokesVector& GetFinalElectronPolarization()
  {
    return fElectronPolarization;
  }
  const G4StokesVector& GetFinalPositronPolarization()
  {
    return fPositronPolarization;
  }

 private:
  G4VPolarizedXS* fCrossSectionCalculator;

  G4StokesVector fBeamPolarization;
  G4StokesVector fTargetPolarization;

  G4StokesVector fPositronPolarization;
  G4StokesVector fElectronPolarization;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
