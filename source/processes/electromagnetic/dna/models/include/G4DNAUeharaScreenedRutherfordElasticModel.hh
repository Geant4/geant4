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

#ifndef G4DNAUeharaScreenedRutherfordElasticModel_h
#define G4DNAUeharaScreenedRutherfordElasticModel_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NistManager.hh"

class G4DNAUeharaScreenedRutherfordElasticModel : public G4VEmModel
{
public:
  G4DNAUeharaScreenedRutherfordElasticModel(const G4ParticleDefinition* p = 0, 
              const G4String& nam = "DNAUeharaScreenedRutherfordElasticModel");

  ~G4DNAUeharaScreenedRutherfordElasticModel() override = default;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double CrossSectionPerVolume(const G4Material* material,
				 const G4ParticleDefinition* p,
				 G4double ekin,
				 G4double emin,
				 G4double emax) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;

  inline void SelectFasterComputation(G4bool input);
  
  //---
  // kept for backward compatibility
  inline void SetKillBelowThreshold (G4double threshold);
  inline G4double GetKillBelowThreshold();
  inline void SelectHighEnergyLimit(G4double threshold);
  //---

  //
  G4DNAUeharaScreenedRutherfordElasticModel&
    operator=(const G4DNAUeharaScreenedRutherfordElasticModel &right) = delete;
  G4DNAUeharaScreenedRutherfordElasticModel(
      const G4DNAUeharaScreenedRutherfordElasticModel&) = delete;

private:

  // -- Cross section
  G4double RutherfordCrossSection(G4double energy, G4double z);
  G4double ScreeningFactor(G4double energy, G4double z);
  
  // -- Final state according to Brenner & Zaider
  G4double BrennerZaiderRandomizeCosTheta(G4double k);
  G4double CalculatePolynomial(G4double k,
                               std::vector<G4double>& vec);
   
  // -- Final state according to Screened Rutherford
  G4double ScreenedRutherfordRandomizeCosTheta(G4double k,
                                               G4double z);

protected:
  G4ParticleChangeForGamma* fParticleChangeForGamma = nullptr;

private:
  G4double iLowEnergyLimit;
  G4double intermediateEnergyLimit;
  G4double iHighEnergyLimit;
  
  // -- Brenner & Zaider
  std::vector<G4double> betaCoeff;
  std::vector<G4double> deltaCoeff;
  std::vector<G4double> gamma035_10Coeff;
  std::vector<G4double> gamma10_100Coeff;
  std::vector<G4double> gamma100_200Coeff;

  // -- Water density table
  const std::vector<G4double>* fpWaterDensity = nullptr;
  
  G4int verboseLevel;
  G4bool isInitialised = false;
  // Selection of computation method
  // We do not recommend "true" usage with the current cumul. proba. settings
  G4bool fasterCode = false;
};
 

inline void G4DNAUeharaScreenedRutherfordElasticModel::
SelectFasterComputation(G4bool input)
{ 
  fasterCode = input; 
}

//---
// kept for backward compatibility

inline void
G4DNAUeharaScreenedRutherfordElasticModel::SelectHighEnergyLimit(
    G4double threshold)
{
  if(threshold > 10. * CLHEP::keV)
  {
    G4Exception (
        "*** WARNING : the G4DNAUeharaScreenedRutherfordElasticModel class is "
        "used above 10 keV !",
        "", JustWarning, "");
  }

  SetHighEnergyLimit(threshold);
}

inline void
G4DNAUeharaScreenedRutherfordElasticModel::SetKillBelowThreshold(G4double)
{
  G4ExceptionDescription errMsg;
  errMsg << "*** WARNING : "
      << "G4DNAUeharaScreenedRutherfordElasticModel::SetKillBelowThreshold"
      << "is deprecated, the kill threshold won't be taken into account";

  G4Exception (
      "G4DNAUeharaScreenedRutherfordElasticModel::SetKillBelowThreshold",
      "DEPRECATED", JustWarning, errMsg);
}

inline G4double
G4DNAUeharaScreenedRutherfordElasticModel::GetKillBelowThreshold()
{
  G4ExceptionDescription errMsg;
  errMsg << "*** WARNING : "
      << "G4DNAUeharaScreenedRutherfordElasticModel::GetKillBelowThreshold"
      << "is deprecated, the returned value is nonsense";

  G4Exception (
      "G4DNAUeharaScreenedRutherfordElasticModel::GetKillBelowThreshold",
      "DEPRECATED", JustWarning, errMsg);

  return -1;
}
//---

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
