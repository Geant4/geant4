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

class G4MollerBhabhaModel : public G4VEmModel
{

public:

  G4MollerBhabhaModel(const G4ParticleDefinition* p = 0);

  ~G4MollerBhabhaModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double HighEnergyLimit(const G4ParticleDefinition* p);

  G4double LowEnergyLimit(const G4ParticleDefinition* p);

  void SetHighEnergyLimit(G4double e) {highKinEnergy = e;};

  void SetLowEnergyLimit(G4double e) {lowKinEnergy = e;};

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple*);

  G4bool IsInCharge(const G4ParticleDefinition*);

  G4double ComputeDEDX(const G4Material*,
                       const G4ParticleDefinition*,
                             G4double kineticEnergy,
                             G4double cutEnergy);

  G4double CrossSection(const G4Material*,
                        const G4ParticleDefinition*,
                              G4double kineticEnergy,
                              G4double cutEnergy,
                              G4double maxEnergy);

  G4DynamicParticle* SampleSecondary(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);

  G4std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);

  G4double MaxSecondaryEnergy(const G4DynamicParticle*);

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                                    G4double kinEnergy);

private:

  void SetParticle(const G4ParticleDefinition* p);

  // hide assignment operator 
  G4MollerBhabhaModel & operator=(const  G4MollerBhabhaModel &right);
  G4MollerBhabhaModel(const  G4MollerBhabhaModel&);

  const G4ParticleDefinition* particle;
  G4ParticleDefinition* theElectron;
  G4double highKinEnergy;
  G4double lowKinEnergy;
  G4double twoln10;
  G4double lowLimit;
  G4bool   isElectron;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MollerBhabhaModel::MaxSecondaryEnergy(
          const G4ParticleDefinition* p,
                G4double kinEnergy) 
{
  G4double tmax = kinEnergy;
  if(isElectron) tmax *= 0.5;
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MollerBhabhaModel::MaxSecondaryEnergy(const G4DynamicParticle* dp)
{
  G4double tmax = dp->GetKineticEnergy();
  if(isElectron) tmax *= 0.5;
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
