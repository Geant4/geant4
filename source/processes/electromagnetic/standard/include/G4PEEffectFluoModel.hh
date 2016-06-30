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
// $Id: G4PEEffectFluoModel.hh 96934 2016-05-18 09:10:41Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PEEffectFluoModel
//
// Author:        Vladimir Ivanchenko on base of G4PEEffectModel
//
// Creation date: 13.06.2010
//
// Modifications:
//
// Class Description:
//
// Implementation of the photo-electric effect with deexcitation
//

// -------------------------------------------------------------------
//

#ifndef G4PEEffectFluoModel_h
#define G4PEEffectFluoModel_h 1

#include "G4VEmModel.hh"
#include <vector>

class G4ParticleChangeForGamma;
class G4VAtomDeexcitation;

class G4PEEffectFluoModel : public G4VEmModel
{

public:

  explicit G4PEEffectFluoModel(const G4String& nam = "PhotoElectric");

  virtual ~G4PEEffectFluoModel();

  virtual 
  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  virtual 
  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
				      G4double kinEnergy,
				      G4double Z,
				      G4double A,
				      G4double, G4double) override;
				      
  virtual G4double CrossSectionPerVolume(const G4Material*,
					 const G4ParticleDefinition*,
					 G4double kineticEnergy,
					 G4double cutEnergy,
					 G4double maxEnergy) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;
private:

  G4PEEffectFluoModel & operator=(const G4PEEffectFluoModel &right) = delete;
  G4PEEffectFluoModel(const G4PEEffectFluoModel&) = delete;

  G4ParticleDefinition*     theGamma;
  G4ParticleDefinition*     theElectron;
  G4ParticleChangeForGamma* fParticleChange;
  G4VAtomDeexcitation*      fAtomDeexcitation;

  G4double                  fminimalEnergy;
  std::vector<G4double>     fSandiaCof;
  std::vector<G4double>     fMatEnergyTh;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
