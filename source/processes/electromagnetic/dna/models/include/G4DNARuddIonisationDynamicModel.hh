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
// Created 11.02.2025 V.Ivanchenko & M. Vologzhin
//                    on base of previous Rudd models
//
// Rudd model of ion ionisation using dynamic mass and charge of an ion
//

#ifndef G4DNARuddIonisationDynamicModel_h
#define G4DNARuddIonisationDynamicModel_h 1

#include "G4VEmModel.hh"
#include "G4VSIntegration.hh"
#include "G4ParticleChangeForGamma.hh"

#include "G4DNAWaterIonisationStructure.hh"
#include <vector>

class G4DNAChemistryManager;
class G4VAtomDeexcitation;
class G4ExtendedPhysicsVector;
class G4Pow;
class G4EmCorrections;

class G4DNARuddIonisationDynamicModel : public G4VEmModel, public G4VSIntegration
{
public:

  explicit G4DNARuddIonisationDynamicModel(const G4ParticleDefinition* p = nullptr,
		           const G4String& nam = "DNARuddIonisationDynamicModel");

  ~G4DNARuddIonisationDynamicModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double ProbabilityDensityFunction(G4double ekin) override;

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

  void StartTracking(G4Track*) override;

  G4DNARuddIonisationDynamicModel & operator=
  (const  G4DNARuddIonisationDynamicModel &right) = delete;
  G4DNARuddIonisationDynamicModel(const G4DNARuddIonisationDynamicModel&) = delete;

private:

  void LoadData();
  
  void SetParticle(const G4ParticleDefinition*);

  G4int SelectShell();

  G4double MaxEnergy();

  G4double SampleElectronEnergy();

  G4double CorrectionFactor();

  G4double S_1s(G4double t,
		G4double energyTransferred,
		G4double slaterEffectiveChg,
		G4double shellNumber);

  G4double S_2s(G4double t,
		G4double energyTransferred,
		G4double slaterEffectiveChg,
		G4double shellNumber);


  G4double S_2p(G4double t,
		G4double energyTransferred,
		G4double slaterEffectiveChg,
		G4double shellNumber);

  G4double Rh(G4double t, 
	      G4double energyTransferred,
	      G4double slaterEffectiveChg,
	      G4double shellNumber);

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma{nullptr};

private:

  // Water density table
  static const std::vector<G4double>* fpWaterDensity;

  // cross section data
  static G4ExtendedPhysicsVector* xsdata_alpha;
  static G4ExtendedPhysicsVector* xsdata_alphap;
  static G4ExtendedPhysicsVector* xsdata_hydrogen;
  static G4ExtendedPhysicsVector* xsdata_helium;
  static G4ExtendedPhysicsVector* xsdata_p;

  // run time data
  G4ExtendedPhysicsVector* xsdata{nullptr};

  const G4ParticleDefinition* fParticle{nullptr};
  const G4Track* fTrack{nullptr};

  G4DNAChemistryManager* fChemistry{nullptr};
  G4EmCorrections* fEmCorrections;
  G4Pow* fGpow;
 
  //deexcitation manager to produce fluo photons and e-
  G4VAtomDeexcitation* fAtomDeexcitation{nullptr};

  // low-energy limit of proton x-section table
  G4double fLowestEnergy{0.0};

  // tracking cut
  G4double fAbsorptionEnergy{0.0};

  G4double fMass{0.0};
  G4double fMassRate{1.0};
  G4double fScaledEnergy{0.0};

  G4double slaterEffectiveCharge[3] = {0.0};
  G4double sCoefficient[3] = {0.0};

  G4double F1{0.0};
  G4double F2{0.0};
  G4double alphaConst{0.0};
  G4double bEnergy{0.0};
  G4double u{0.0};
  G4double v{0.0};
  G4double wc{0.0};

  G4int fSelectedShell{0};
  G4int verbose{0};
  std::size_t idx{0};

  G4bool isFirst{false};
  G4bool isInitialised{false};
  G4bool isIon{false};
  G4bool isHelium{false};
  G4bool statCode{false};
  G4bool useDNAWaterStructure{true};

  // energy levels of water molecule  
  G4DNAWaterIonisationStructure waterStructure;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
