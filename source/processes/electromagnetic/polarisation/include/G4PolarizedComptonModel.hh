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
// File name:     G4PolarizedComptonModel
//
// Author:        Andreas Schaelicke
//
// Class Description:
//   Implementation of polarized gamma Compton scattering on free electron

#ifndef G4PolarizedComptonModel_h
#define G4PolarizedComptonModel_h 1

#include "globals.hh"
#include "G4KleinNishinaCompton.hh"
#include "G4StokesVector.hh"
#include "G4ThreeVector.hh"

class G4DynamicParticle;
class G4MaterialCutsCouple;
class G4ParticleDefinition;
class G4ParticleChangeForGamma;
class G4VPolarizedXS;

class G4PolarizedComptonModel : public G4KleinNishinaCompton
{
 public:
  explicit G4PolarizedComptonModel(const G4ParticleDefinition* p = nullptr,
                                   const G4String& nam = "Polarized-Compton");

  virtual ~G4PolarizedComptonModel() override;
  
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                              G4double kinEnergy, G4double Z,
                                              G4double A, G4double cut,
                                              G4double emax) override;
  G4double ComputeAsymmetryPerAtom(G4double gammaEnergy, G4double Z);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*, G4double tmin,
                                 G4double maxEnergy) override;

  // polarized routines
  inline void SetTargetPolarization(const G4ThreeVector& pTarget);
  inline void SetBeamPolarization(const G4ThreeVector& pBeam);
  inline const G4ThreeVector& GetTargetPolarization() const;
  inline const G4ThreeVector& GetBeamPolarization() const;
  inline const G4ThreeVector& GetFinalGammaPolarization() const;
  inline const G4ThreeVector& GetFinalElectronPolarization() const;

  // hide assignment operator
  G4PolarizedComptonModel& operator=(const G4PolarizedComptonModel& right) =
    delete;
  G4PolarizedComptonModel(const G4PolarizedComptonModel&) = delete;

 private:
  void PrintWarning(const G4DynamicParticle*, G4int, G4double grej,
                    G4double onecos, G4double phi, const G4String) const;

  static constexpr G4int fLoopLim = 10000;

  G4VPolarizedXS* fCrossSectionCalculator;
  // incoming
  G4StokesVector fBeamPolarization;    // photon
  G4StokesVector fTargetPolarization;  // electron
  // outgoing
  G4StokesVector fFinalGammaPolarization;
  G4StokesVector finalElectronPolarization;

  G4int fVerboseLevel;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4PolarizedComptonModel::SetTargetPolarization(
  const G4ThreeVector& pTarget)
{
  fTargetPolarization = G4StokesVector(pTarget);
}
inline void G4PolarizedComptonModel::SetBeamPolarization(
  const G4ThreeVector& pBeam)
{
  fBeamPolarization = G4StokesVector(pBeam);
}
inline const G4ThreeVector& G4PolarizedComptonModel::GetTargetPolarization()
  const
{
  return fTargetPolarization;
}
inline const G4ThreeVector& G4PolarizedComptonModel::GetBeamPolarization() const
{
  return fBeamPolarization;
}
inline const G4ThreeVector& G4PolarizedComptonModel::GetFinalGammaPolarization()
  const
{
  return fFinalGammaPolarization;
}
inline const G4ThreeVector&
G4PolarizedComptonModel::GetFinalElectronPolarization() const
{
  return finalElectronPolarization;
}

#endif
