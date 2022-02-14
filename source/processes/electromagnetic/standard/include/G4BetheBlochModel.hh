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
// File name:     G4BetheBlochModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 03.01.2002
//
// Modifications:
//
// Modified by Michel Maire and Vladimir Ivanchenko
//
// Class Description:
//
// Implementation of Bethe-Bloch model of energy loss and
// delta-electron production by heavy charged particles

// -------------------------------------------------------------------
//

#ifndef G4BetheBlochModel_h
#define G4BetheBlochModel_h 1

#include "G4VEmModel.hh"

class G4EmCorrections;
class G4ParticleChangeForLoss;
class G4ICRU90StoppingData;
class G4NistManager;

class G4BetheBlochModel : public G4VEmModel
{

public:

  explicit G4BetheBlochModel(const G4ParticleDefinition* p = nullptr,
			     const G4String& nam = "BetheBloch");

  ~G4BetheBlochModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double MinEnergyCut(const G4ParticleDefinition*,
			const G4MaterialCutsCouple* couple) override;

  virtual G4double ComputeCrossSectionPerElectron(
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

  G4double GetChargeSquareRatio(const G4ParticleDefinition* p,
				const G4Material* mat,
				G4double kineticEnergy) override;

  G4double GetParticleCharge(const G4ParticleDefinition* p,
			     const G4Material* mat,
			     G4double kineticEnergy) override;

  void CorrectionsAlongStep(const G4MaterialCutsCouple* couple,
			    const G4DynamicParticle* dp,
			    const G4double& length,
			    G4double& eloss) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;

  // hide assignment operator
  G4BetheBlochModel & operator=(const G4BetheBlochModel &right) = delete;
  G4BetheBlochModel(const G4BetheBlochModel&) = delete;

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
			      G4double kinEnergy) override;

  inline G4double GetChargeSquareRatio() const;

  inline void SetChargeSquareRatio(G4double val);

private:

  void SetupParameters(const G4ParticleDefinition* p);

  const G4ParticleDefinition* particle = nullptr;
  const G4ParticleDefinition* theElectron;
  G4EmCorrections*            corr;
  G4ParticleChangeForLoss*    fParticleChange = nullptr;
  G4NistManager*              nist;
  G4ICRU90StoppingData*       fICRU90 = nullptr;
  const G4Material*           currentMaterial = nullptr;
  const G4Material*           baseMaterial = nullptr;

  G4double mass = 0.0;
  G4double tlimit = DBL_MAX;
  G4double spin = 0.0;
  G4double magMoment2 = 0.0;
  G4double chargeSquare = 1.0;
  G4double ratio = 1.0;
  G4double formfact = 0.0;
  G4double twoln10;
  G4double fAlphaTlimit;
  G4double fProtonTlimit;

  G4int    iICRU90 = -1;
  G4bool   isIon = false;
  G4bool   isAlpha = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4BetheBlochModel::GetChargeSquareRatio() const
{
  return chargeSquare;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4BetheBlochModel::SetChargeSquareRatio(G4double val)
{
  chargeSquare = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
