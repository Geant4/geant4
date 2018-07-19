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
// $Id: G4BraggModel.hh 96934 2016-05-18 09:10:41Z gcosmo $
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

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VEmModel.hh"
#include "G4PSTARStopping.hh"

class G4ParticleChangeForLoss;
class G4EmCorrections;

class G4BraggModel : public G4VEmModel
{

public:

  explicit G4BraggModel(const G4ParticleDefinition* p = nullptr,
			const G4String& nam = "Bragg");

  virtual ~G4BraggModel();

  virtual void Initialise(const G4ParticleDefinition*, 
			  const G4DataVector&) override;

  virtual G4double ComputeCrossSectionPerElectron(
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double cutEnergy,
				 G4double maxEnergy);
				 
  virtual G4double ComputeCrossSectionPerAtom(
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double Z, G4double A,
				 G4double cutEnergy,
				 G4double maxEnergy) override;
				 				 
  virtual G4double CrossSectionPerVolume(const G4Material*,
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double cutEnergy,
				 G4double maxEnergy) override;
				 
  virtual G4double ComputeDEDXPerVolume(const G4Material*,
				const G4ParticleDefinition*,
				G4double kineticEnergy,
				G4double cutEnergy) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;

  // Compute ion charge 
  virtual G4double GetChargeSquareRatio(const G4ParticleDefinition*,
					const G4Material*,
					G4double kineticEnergy) override;

  virtual G4double GetParticleCharge(const G4ParticleDefinition* p,
				     const G4Material* mat,
				     G4double kineticEnergy) override;

protected:

  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
				      G4double kinEnergy) final;

  inline G4double GetChargeSquareRatio() const;

  inline void SetChargeSquareRatio(G4double val);

private:

  inline void SetParticle(const G4ParticleDefinition* p);

  G4bool HasMaterial(const G4Material* material);

  G4double StoppingPower(const G4Material* material,
                               G4double kineticEnergy);

  G4double ElectronicStoppingPower(G4double z,
                                   G4double kineticEnergy) const;

  G4double DEDX(const G4Material* material, G4double kineticEnergy);

  G4bool MolecIsInZiegler1988(const G4Material* material);

  G4double ChemicalFactor(G4double kineticEnergy, G4double eloss125) const;

  // hide assignment operator
  G4BraggModel & operator=(const  G4BraggModel &right) = delete;
  G4BraggModel(const  G4BraggModel&) = delete;

  G4EmCorrections*            corr;

  const G4ParticleDefinition* particle;
  G4ParticleDefinition*       theElectron;
  G4ParticleChangeForLoss*    fParticleChange;
  static G4PSTARStopping*     fPSTAR;

  const G4Material*           currentMaterial;

  G4double mass;
  G4double spin;
  G4double chargeSquare;
  G4double massRate;
  G4double ratio;
  G4double lowestKinEnergy;
  G4double protonMassAMU;
  G4double theZieglerFactor;
  G4double expStopPower125;    // Experimental Stopping power at 125keV

  G4int    iMolecula;          // index in the molecula's table
  G4int    iPSTAR;             // index in PSTAR
  G4bool   isIon;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4BraggModel::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  mass = particle->GetPDGMass();
  spin = particle->GetPDGSpin();
  G4double q = particle->GetPDGCharge()/CLHEP::eplus;
  chargeSquare = q*q;
  massRate     = mass/CLHEP::proton_mass_c2;
  ratio = CLHEP::electron_mass_c2/mass;
}

inline G4double G4BraggModel::GetChargeSquareRatio() const
{
  return chargeSquare;
}

inline void G4BraggModel::SetChargeSquareRatio(G4double val)
{
  chargeSquare = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
