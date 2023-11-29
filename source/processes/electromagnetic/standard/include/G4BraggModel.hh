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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BraggModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications:
// 23-12-02 V.Ivanchenko change interface in order to moveto cut per region
// 24-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 12-11-03 Fix for GenericIons (V.Ivanchenko)
// 11-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 15-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 25-04-06 Added stopping data from PSTAR (V.Ivanchenko)
// 12-08-08 Added methods GetParticleCharge, GetChargeSquareRatio, 
//          CorrectionsAlongStep needed for ions(V.Ivanchenko)

//
// Class Description:
//
// Implementation of energy loss and delta-electron production
// by heavy slow charged particles using ICRU'49 and NIST evaluated data 
// for protons

// -------------------------------------------------------------------
//

#ifndef G4BraggModel_h
#define G4BraggModel_h 1

#include "G4VEmModel.hh"
#include "G4PSTARStopping.hh"

class G4ParticleChangeForLoss;
class G4EmCorrections;
class G4ICRU90StoppingData;
class G4PSTARStopping;

class G4BraggModel : public G4VEmModel
{

public:

  explicit G4BraggModel(const G4ParticleDefinition* p = nullptr,
			const G4String& nam = "Bragg");

  ~G4BraggModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double MinEnergyCut(const G4ParticleDefinition*,
			const G4MaterialCutsCouple* couple) override;

  G4double ComputeCrossSectionPerElectron(
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double cutEnergy,
				 G4double maxEnergy);
				 
  G4double ComputeCrossSectionPerAtom(
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double Z, G4double A,
				 G4double cutEnergy,
				 G4double maxEnergy) override;
				 				 
  G4double CrossSectionPerVolume(const G4Material*,
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double cutEnergy,
				 G4double maxEnergy) override;
				 
  G4double ComputeDEDXPerVolume(const G4Material*,
				const G4ParticleDefinition*,
				G4double kineticEnergy,
				G4double cutEnergy) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;

  // Compute ion charge 
  G4double GetChargeSquareRatio(const G4ParticleDefinition*,
				const G4Material*,
				G4double kineticEnergy) override;

  G4double GetParticleCharge(const G4ParticleDefinition* p,
			     const G4Material* mat,
			     G4double kineticEnergy) override;

  // hide assignment operator
  G4BraggModel & operator=(const  G4BraggModel &right) = delete;
  G4BraggModel(const  G4BraggModel&) = delete;

protected:

  void SetParticle(const G4ParticleDefinition* p);

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
			      G4double kinEnergy) final;

  inline void SetChargeSquareRatio(G4double val);

private:

  void HasMaterial(const G4Material* material);

  G4double StoppingPower(const G4Material* material,
			 G4double kineticEnergy);

  G4double ElectronicStoppingPower(G4double z,
                                   G4double kineticEnergy) const;

  G4double DEDX(const G4Material* material, G4double kineticEnergy);

  G4bool MolecIsInZiegler1988(const G4Material* material);

  G4double ChemicalFactor(G4double kineticEnergy, G4double eloss125) const;

protected:

  const G4ParticleDefinition* particle = nullptr;
  G4ParticleDefinition* theElectron = nullptr;
  G4ParticleChangeForLoss* fParticleChange = nullptr;

  const G4Material* currentMaterial = nullptr;
  const G4Material* baseMaterial = nullptr;

  G4EmCorrections* corr = nullptr;

  static G4ICRU90StoppingData* fICRU90;
  static G4PSTARStopping* fPSTAR;

  G4double mass = 0.0;
  G4double spin = 0.0;
  G4double chargeSquare = 1.0;
  G4double massRate = 1.0;
  G4double ratio = 1.0;
  G4double protonMassAMU = 1.007276;
  G4double lowestKinEnergy;
  G4double theZieglerFactor;
  G4double expStopPower125;  // Experimental Stopping power at 125keV

  G4int iMolecula = -1;   // index in the molecula's table
  G4int iPSTAR = -1;      // index in NIST PSTAR
  G4int iICRU90 = -1;

private:

  G4bool isIon = false;
  G4bool isFirst = false;
};

inline void G4BraggModel::SetChargeSquareRatio(G4double val)
{
  chargeSquare = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
