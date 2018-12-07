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
// File name:     G4eplusTo2GammaOKVIModel
//
// Author:        Vladimir Ivanchenko and Omrane Kadri 
//
// Creation date: 29.03.2018
//
//
// Class Description:
//
// Implementation of e+ annihilation into 2 gamma

// -------------------------------------------------------------------
//

#ifndef G4eplusTo2GammaOKVIModel_h
#define G4eplusTo2GammaOKVIModel_h 1

#include "G4VEmModel.hh"

class G4eplusTo3GammaOKVIModel;
class G4ParticleChangeForGamma;
class G4PhysicsVector;
class G4DataVector;

class G4eplusTo2GammaOKVIModel : public G4VEmModel
{

public:

  explicit G4eplusTo2GammaOKVIModel(const G4ParticleDefinition* p = nullptr,
                                    const G4String& nam = "eplus2ggOKVI");

  virtual ~G4eplusTo2GammaOKVIModel();

  virtual void Initialise(const G4ParticleDefinition*, 
			  const G4DataVector&) override;

  G4double ComputeCrossSectionPerElectron(G4double kinEnergy); 
  
  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A = 0., 
                                      G4double cutEnergy = 0.,
                                      G4double maxEnergy = DBL_MAX) final;

  virtual G4double CrossSectionPerVolume(const G4Material*,
					 const G4ParticleDefinition*,
					 G4double kineticEnergy,
					 G4double cutEnergy = 0.0,
					 G4double maxEnergy = DBL_MAX) final;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin = 0.0,
				 G4double maxEnergy = DBL_MAX) final;

  inline void SetDelta(G4double val) { if(val > 0.0) { fDelta = val; } };

private:

  // hide assignment operator
  G4eplusTo2GammaOKVIModel & operator=
  (const  G4eplusTo2GammaOKVIModel &right) = delete;
  G4eplusTo2GammaOKVIModel(const  G4eplusTo2GammaOKVIModel&) = delete;

  G4double fDelta;
  G4double fGammaTh;
  const G4ParticleDefinition* theGamma;
  G4ParticleChangeForGamma* fParticleChange;
  G4eplusTo3GammaOKVIModel* f3GModel;

  const G4DataVector*     fCuts;
 
  static G4PhysicsVector* fCrossSection;
  static G4PhysicsVector* fCrossSection3G;
  static G4PhysicsVector* f3GProbability;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
