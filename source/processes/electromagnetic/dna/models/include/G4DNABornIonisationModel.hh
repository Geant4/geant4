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
// Created 25.03.2025 V.Ivanchenko 
// on base of the G4DNABornIonisationModel1 of S.Incerti & M.Karamitros
//
// Simulation of ionisation for electrons and protons
//

#ifndef G4DNABornIonisationModel_h
#define G4DNABornIonisationModel_h 1

#include "G4VEmModel.hh"
#include "G4VSIntegration.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4DNAWaterIonisationStructure.hh"

class G4DNAChemistryManager;
class G4VAtomDeexcitation;
class G4DNACrossSectionDataSet;
class G4DNASamplingTable;

class G4DNABornIonisationModel : public G4VEmModel, public G4VSIntegration
{
public:

  G4DNABornIonisationModel(const G4ParticleDefinition* p = nullptr,
		           const G4String& nam = "DNABornIonisationModel");

  ~G4DNABornIonisationModel() override;
   
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
			 G4double tmin, G4double maxEnergy) override;

  void StartTracking(G4Track*) override;
  
  void SelectFasterComputation(G4bool input) { fasterCode = input; }; 

  void SelectStationary(G4bool input) { statCode = input; }; 

  void SelectSPScaling(G4bool input) { spScaling = input; };

  G4DNABornIonisationModel & operator=(const  G4DNABornIonisationModel &right) = delete;
  G4DNABornIonisationModel(const  G4DNABornIonisationModel&) = delete;

private:

  void LoadData();

  G4int SelectShell();

  G4double SampleCumulative();

  G4double SampleDifferential();

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:

  // Water density table
  static const std::vector<G4double>* fpWaterDensity;

  // data
  static G4DNACrossSectionDataSet* xsdata_e;
  static G4DNACrossSectionDataSet* xsdata_p;
  G4DNACrossSectionDataSet* xsdata{nullptr};

  // sampling data
  static G4DNASamplingTable* sampling_e;
  static G4DNASamplingTable* sampling_p;
  G4DNASamplingTable* sampling;
    
  const G4ParticleDefinition* fParticle{nullptr};
  const G4Track* fTrack{nullptr};

  G4DNAChemistryManager* fChemistry{nullptr};

  // Deexcitation manager to produce fluo photons and e-
  G4VAtomDeexcitation* fAtomDeexcitation;

  // limits of x-section table
  G4double fLowEnergy{0.0};
  G4double fHighEnergy{0.0};
  G4double fpLimitEnergy{0.0};
  G4double feLimitEnergy{0.0};

  // tracking cut
  G4double fAbsorptionEnergy{0.0};

  G4double fMass{0.0};
  G4double fPrimaryEnergy{0.0};
  G4double fMaxEnergy{0.0};
  G4double fTemp[5] = {0.0};

  G4int fSelectedShell{0};
  G4int verbose{0};

  G4bool isFirst{false};
  G4bool isInitialised{false};
  G4bool isElectron{false};
  G4bool fasterCode{false};
  G4bool statCode{false};
  G4bool spScaling{true};

  // Final state  
  G4DNAWaterIonisationStructure waterStructure;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
