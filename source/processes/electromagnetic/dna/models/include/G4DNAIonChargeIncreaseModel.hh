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
// Author V.Ivanchenko 15.04.2024
//
// General ion model for simulation of charge increase.
// Concrete model of ion ionisation is selected on fly.
//

#ifndef G4DNAIonChargeIncreaseModel_h
#define G4DNAIonChargeIncreaseModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"

class G4DNAIonChargeIncreaseModel : public G4VEmModel
{
public:

  G4DNAIonChargeIncreaseModel(const G4ParticleDefinition* p = nullptr, 
		              const G4String& nam = "DNAIonChargeIncrease");

  ~G4DNAIonChargeIncreaseModel() override = default;

  G4DNAIonChargeIncreaseModel & operator=
  (const  G4DNAIonChargeIncreaseModel &right) = delete;
  G4DNAIonChargeIncreaseModel(const G4DNAIonChargeIncreaseModel&) = delete;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double CrossSectionPerVolume(const G4Material* material,
				 const G4ParticleDefinition* p,
				 G4double ekin, G4double,
				 G4double) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double, G4double) override;

  void StartTracking(G4Track*) override;

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma{nullptr};

private:

  const G4DynamicParticle* fDynParticle{nullptr};
  G4VEmModel* fCurrentModel{nullptr};

  // list of concrete ion ionisation models, put extra here
  G4VEmModel* fDummy;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
