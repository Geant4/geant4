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
// $Id: G4MollerBhabhaModel.hh 96934 2016-05-18 09:10:41Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4MollerBhabhaModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 07.01.2002
//
// Modifications:
//
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 06-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 14-01-07 promote SetParticle() from private to protected (mma)

//
// Class Description:
//
// Implementation of energy loss process and delta-electrons production
// of electrons and positrons

// -------------------------------------------------------------------
//

#ifndef G4MollerBhabhaModel_h
#define G4MollerBhabhaModel_h 1

#include "G4VEmModel.hh"

class G4ParticleChangeForLoss;

class G4MollerBhabhaModel : public G4VEmModel
{

public:

  explicit G4MollerBhabhaModel(const G4ParticleDefinition* p = nullptr, 
			       const G4String& nam = "MollerBhabha");

  virtual ~G4MollerBhabhaModel();

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

protected:

  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
				      G4double kinEnergy) final;
			      
  inline void SetParticle(const G4ParticleDefinition* p);

  const G4ParticleDefinition* particle;
  G4ParticleDefinition*       theElectron;
  G4ParticleChangeForLoss*    fParticleChange;

  G4bool   isElectron;
  G4double twoln10;
  G4double lowLimit;
  
private:

  // hide assignment operator 
  G4MollerBhabhaModel & operator=(const  G4MollerBhabhaModel &right) = delete;
  G4MollerBhabhaModel(const  G4MollerBhabhaModel&) = delete;

  G4bool   isInitialised;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4MollerBhabhaModel::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  if(p != theElectron) { isElectron = false; }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
