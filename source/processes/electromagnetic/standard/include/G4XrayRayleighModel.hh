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
// author: Vladimir Grichine
//
// 14.10.12 V. Grichine, fFormFactor was added as the class member
// 25.05.2011 first implementation
//
// X ray Rayleigh scattering model based on simplified form-factors 
// and angular distribution
//

#ifndef G4XrayRayleighModel_h
#define G4XrayRayleighModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Gamma.hh"

class G4XrayRayleighModel : public G4VEmModel
{

public:

  explicit G4XrayRayleighModel(const G4ParticleDefinition* p = nullptr,
			       const G4String& nam = "XrayRayleigh");

  ~G4XrayRayleighModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A=0, 
                                      G4double cut=0,
                                      G4double emax=DBL_MAX) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;

  G4XrayRayleighModel & operator=(const  G4XrayRayleighModel &right) = delete;
  G4XrayRayleighModel(const  G4XrayRayleighModel&) = delete;

private:

  G4ParticleChangeForGamma* fParticleChange;

  G4double lowEnergyLimit;  
  G4double highEnergyLimit; 
  G4double fFormFactor; 

  G4int verboseLevel;
  G4bool isInitialised;

  static const G4double fCofA;
  static const G4double fCofR;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
