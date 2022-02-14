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
// File name:     G4PolarizedAnnihilationModel
//
// Author:        Andreas Schaelicke and Pavel Starovoitov
//
// Class Description:
//   Implementation of polarized gamma Annihilation scattering on free electron
//
// -------------------------------------------------------------------

#ifndef G4PolarizedAnnihilationModel_h
#define G4PolarizedAnnihilationModel_h 1

#include "globals.hh"
#include "G4eeToTwoGammaModel.hh"
#include "G4ThreeVector.hh"
#include "G4StokesVector.hh"

class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4ParticleChangeForGamma;
class G4ParticleDefinition;
class G4PolarizedAnnihilationXS;

class G4PolarizedAnnihilationModel : public G4eeToTwoGammaModel
{
 public:
  explicit G4PolarizedAnnihilationModel(
    const G4ParticleDefinition* p = nullptr,
    const G4String& nam           = "Polarized-Annihilation");

  virtual ~G4PolarizedAnnihilationModel() override;

  virtual void Initialise(const G4ParticleDefinition*,
                          const G4DataVector&) override;

  virtual G4double ComputeCrossSectionPerElectron(G4double kinEnergy) override;

  void ComputeAsymmetriesPerElectron(G4double gammaEnergy, G4double& valueX,
                                     G4double& valueA, G4double& valueT);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*, G4double tmin,
                                 G4double maxEnergy) final;

  // polarized routines
  inline void SetTargetPolarization(const G4ThreeVector& pTarget);
  inline void SetBeamPolarization(const G4ThreeVector& pBeam);
  inline const G4ThreeVector& GetTargetPolarization() const;
  inline const G4ThreeVector& GetBeamPolarization() const;
  inline const G4ThreeVector& GetFinalGamma1Polarization() const;
  inline const G4ThreeVector& GetFinalGamma2Polarization() const;

  G4PolarizedAnnihilationModel& operator                            =(
    const G4PolarizedAnnihilationModel& right) = delete;
  G4PolarizedAnnihilationModel(const G4PolarizedAnnihilationModel&) = delete;

 private:
  G4PolarizedAnnihilationXS* fCrossSectionCalculator;
  G4ParticleChangeForGamma* fParticleChange;

  // incoming
  G4StokesVector fBeamPolarization;    // positron
  G4StokesVector fTargetPolarization;  // electron
  // outgoing
  G4StokesVector fFinalGamma1Polarization;
  G4StokesVector fFinalGamma2Polarization;

  G4int fVerboseLevel;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4PolarizedAnnihilationModel::SetTargetPolarization(
  const G4ThreeVector& pTarget)
{
  fTargetPolarization = G4StokesVector(pTarget);
}
inline void G4PolarizedAnnihilationModel::SetBeamPolarization(
  const G4ThreeVector& pBeam)
{
  fBeamPolarization = G4StokesVector(pBeam);
}
inline const G4ThreeVector&
G4PolarizedAnnihilationModel::GetTargetPolarization() const
{
  return fTargetPolarization;
}
inline const G4ThreeVector& G4PolarizedAnnihilationModel::GetBeamPolarization()
  const
{
  return fBeamPolarization;
}
inline const G4ThreeVector&
G4PolarizedAnnihilationModel::GetFinalGamma1Polarization() const
{
  return fFinalGamma1Polarization;
}
inline const G4ThreeVector&
G4PolarizedAnnihilationModel::GetFinalGamma2Polarization() const
{
  return fFinalGamma2Polarization;
}

#endif
