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
// Created by Z. Francis, S. Incerti 2010 
//
// Modified for inverse rudd function sampling 26.10.2010
// Rewitten by V.Ivanchenko 21.05.2023
//

#ifndef G4DNARuddIonisationExtendedModel_h
#define G4DNARuddIonisationExtendedModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"

#include "G4EmCorrections.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4DNACrossSectionDataSet.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4LogLogInterpolation.hh"

#include "G4DNAWaterIonisationStructure.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4NistManager.hh"
#include <vector>

class G4DNARuddIonisationExtendedModel : public G4VEmModel
{
public:

  explicit G4DNARuddIonisationExtendedModel(const G4ParticleDefinition* p = nullptr,
		           const G4String& nam = "DNARuddIonisationExtendedModel");

  ~G4DNARuddIonisationExtendedModel() override;

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

  void SelectStationary(G4bool val) { statCode = val; }; 

  // method for unit tests
  G4double ComputeProbabilityFunction(const G4ParticleDefinition*, G4double kine,
                                      G4double deltae, G4int shell);

  G4DNARuddIonisationExtendedModel & operator=
  (const  G4DNARuddIonisationExtendedModel &right) = delete;
  G4DNARuddIonisationExtendedModel(const G4DNARuddIonisationExtendedModel&) = delete;

private:

  void SetParticle(const G4ParticleDefinition*);

  G4int SelectShell(G4double energy);

  G4double SampleElectronEnergy(G4double kine, G4double bindingEnergy, G4int shell);

  G4double ProbabilityFunction(G4double kine, G4double deltae,
                               G4double bindingEnergy, G4int shell);

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

  G4double CorrectionFactor(G4double kine, G4int shell);

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma{nullptr};

private:

  // idx = 0 - hydrogen
  // idx = 1 - proton
  // idx = 2%26 - ions idx=Z
  // idx = -1 - alpha+
  // idx = -1 - helium
  static const G4int RUDDZMAX = 27;
  static G4DNACrossSectionDataSet* xsdata[RUDDZMAX];
  static G4DNACrossSectionDataSet* xsalphaplus;
  static G4DNACrossSectionDataSet* xshelium;

  // Water density table
  static const std::vector<G4double>* fpWaterDensity;

  G4DNACrossSectionDataSet* xscurrent{nullptr};
  const G4ParticleDefinition* fParticle{nullptr};
  G4EmCorrections* fEmCorrections;
  G4Pow* fGpow;
 
  //deexcitation manager to produce fluo photons and e-
  G4VAtomDeexcitation* fAtomDeexcitation{nullptr};

  // tracking cut and low-energy limit of proton x-section table
  G4double fLowestEnergy;
  // scaled low-energy limit of ion x-section table
  G4double fLimitEnergy;

  G4double fMass{0.0};
  G4double fAmass{0.0};
  G4double fMassRate{1.0};
  G4double fElow{0.0};

  G4double slaterEffectiveCharge[3] = {0.0};
  G4double sCoefficient[3] = {0.0};
  G4double fTemp[5] = {0.0};

  G4int idx{-1};
  G4int verbose{0};

  G4bool isIon{false};
  G4bool isFirst{false};
  G4bool isHelium{false};
  G4bool statCode{false};
  G4bool useDNAWaterStructure{true};

  // energy levels of water molecule  
  G4DNAWaterIonisationStructure waterStructure;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
