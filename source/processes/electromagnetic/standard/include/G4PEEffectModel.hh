//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4PEEffectModel.hh,v 1.3 2005/05/12 11:06:43 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PEEffectModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 21.04.2005
//
// Modifications:
//
// Class Description:
//
// Implementation of the photo-electric effect
//

// -------------------------------------------------------------------
//

#ifndef G4PEEffectModel_h
#define G4PEEffectModel_h 1

#include "G4VEmModel.hh"
#include "G4PhysicsTable.hh"

class G4ParticleChangeForGamma;

class G4PEEffectModel : public G4VEmModel
{

public:

  G4PEEffectModel(const G4ParticleDefinition* p = 0,
		  const G4String& nam = "PhotoElectric");

  virtual ~G4PEEffectModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy,
                                      G4double Z,
                                      G4double A,
                                      G4double cut,
                                      G4double emax);

  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);
protected:

  virtual G4double ElecCosThetaDistribution(G4double ElecKineEnergy);

private:

  G4ParticleDefinition*     theGamma;
  G4ParticleDefinition*     theElectron;
  G4ParticleChangeForGamma* fParticleChange;

  G4double                  fminimalEnergy;
  G4bool                    isInitialized;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

inline G4double G4PEEffectModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double energy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
 G4double* SandiaCof = G4SandiaTable::GetSandiaCofPerAtom((G4int)Z, energy);

 G4double energy2 = energy*energy, energy3 = energy*energy2,
          energy4 = energy2*energy2;

 return SandiaCof[0]/energy  + SandiaCof[1]/energy2 +
        SandiaCof[2]/energy3 + SandiaCof[3]/energy4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
