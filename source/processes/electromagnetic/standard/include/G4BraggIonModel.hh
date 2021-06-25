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
// File name:     G4BraggIonModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.10.2004
//
// Modifications:
// 09-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 11-05-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 15-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 25-04-06 Add stopping data from ASTAR (V.Ivanchenko)
// 12-08-08 Added methods GetParticleCharge, GetChargeSquareRatio, 
//          CorrectionsAlongStep needed for ions(V.Ivanchenko)

//
// Class Description:
//
// Implementation of energy loss and delta-electron production
// by heavy slow charged particles using ICRU'49 and NIST evaluated data 
// for He4 ions

// -------------------------------------------------------------------
//

#ifndef G4BraggIonModel_h
#define G4BraggIonModel_h 1

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VEmModel.hh"
#include "G4ASTARStopping.hh"

class G4ParticleChangeForLoss;
class G4EmCorrections;
class G4ICRU90StoppingData;

class G4BraggIonModel : public G4VEmModel
{

public:

  explicit G4BraggIonModel(const G4ParticleDefinition* p = nullptr,
			   const G4String& nam = "BraggIon");

  ~G4BraggIonModel() override;

  void Initialise(const G4ParticleDefinition*, 
		  const G4DataVector&) override;

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

  // add correction to energy loss and ompute non-ionizing energy loss
  void CorrectionsAlongStep(const G4MaterialCutsCouple*,
			    const G4DynamicParticle*,
			    const G4double& length,
			    G4double& eloss) override;

  // hide assignment operator
  G4BraggIonModel & operator=(const  G4BraggIonModel &right) = delete;
  G4BraggIonModel(const  G4BraggIonModel&) = delete;

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
			      G4double kinEnergy) final;

private:

  void SetParticle(const G4ParticleDefinition* p);

  G4double HeEffChargeSquare(G4double z, G4double kinEnergyInMeV) const;

  void HasMaterial(const G4Material* material);

  G4double StoppingPower(const G4Material* material,
                               G4double kineticEnergy);

  G4double ElectronicStoppingPower(G4double z,
                                   G4double kineticEnergy) const;

  G4double DEDX(const G4Material* material, G4double kineticEnergy);

  G4EmCorrections*            corr = nullptr;
  const G4ParticleDefinition* particle = nullptr;
  G4ParticleDefinition*       theElectron = nullptr;
  G4ParticleChangeForLoss*    fParticleChange = nullptr;
  G4ICRU90StoppingData*       fICRU90 = nullptr;

  const G4Material*           currentMaterial = nullptr;
  const G4Material*           baseMaterial = nullptr;

  static G4ASTARStopping*      fASTAR;

  G4double mass;
  G4double spin;
  G4double chargeSquare;
  G4double massRate;
  G4double ratio;
  G4double lowestKinEnergy;
  G4double HeMass;
  G4double massFactor;
  G4double corrFactor = 1.0;
  G4double rateMassHe2p;
  G4double theZieglerFactor;

  G4int    iMolecula = -1; // index in the molecula's table
  G4int    iASTAR = -1;    // index in ASTAR
  G4int    iICRU90 = -1;
  G4bool   isIon = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4BraggIonModel::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  mass = particle->GetPDGMass();
  spin = particle->GetPDGSpin();
  G4double q   = particle->GetPDGCharge()/CLHEP::eplus;
  chargeSquare = q*q;
  massRate     = mass/CLHEP::proton_mass_c2;
  ratio        = CLHEP::electron_mass_c2/mass;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
