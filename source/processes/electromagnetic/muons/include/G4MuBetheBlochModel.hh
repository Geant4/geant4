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
// File name:     G4MuBetheBlochModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 09.08.2002
//
// Modifications:
//
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add Nama (V.Ivanchenko)
// 10-02-04 Calculation of radiative corrections using R.Kokoulin model (V.I)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 12-04-05 Add usage of G4EmCorrections (V.Ivanchenko)
// 13-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
//

//
// Class Description:
//
// Implementation of Bethe-Bloch model of energy loss and
// delta-electron production by heavy charged particles

// -------------------------------------------------------------------
//

#ifndef G4MuBetheBlochModel_h
#define G4MuBetheBlochModel_h 1

#include "G4VEmModel.hh"

class G4ParticleChangeForLoss;
class G4EmCorrections;

class G4MuBetheBlochModel : public G4VEmModel
{

public:

  explicit G4MuBetheBlochModel(const G4ParticleDefinition* p = nullptr,
                               const G4String& nam = "MuBetheBloch");

  ~G4MuBetheBlochModel() = default;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double MinEnergyCut(const G4ParticleDefinition*,
			const G4MaterialCutsCouple*) override;
			
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

  // hide assignment operator
  G4MuBetheBlochModel & operator=(const  G4MuBetheBlochModel &right) = delete;
  G4MuBetheBlochModel(const  G4MuBetheBlochModel&) = delete;

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                              G4double kinEnergy) override;

private:

  void SetParticle(const G4ParticleDefinition* p);

  const G4ParticleDefinition* particle = nullptr;
  G4ParticleDefinition*       theElectron = nullptr;
  G4ParticleChangeForLoss*    fParticleChange = nullptr;
  G4EmCorrections*            corr = nullptr;

  G4double limitKinEnergy;
  G4double logLimitKinEnergy;
  G4double mass = 1.0;
  G4double massSquare = 1.0;
  G4double ratio = 1.0;
  G4double twoln10;
  G4double alphaprime;
  static G4double xgi[8],wgi[8];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif
