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

class G4MuBetheBlochModel : public G4VEmModel
{

public:

  G4MuBetheBlochModel(const G4ParticleDefinition* p = 0, const G4String& nam = "MuBetheBloch");

  ~G4MuBetheBlochModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double HighEnergyLimit(const G4ParticleDefinition* p);

  G4double LowEnergyLimit(const G4ParticleDefinition* p);

  void SetHighEnergyLimit(G4double e) {highKinEnergy = e;};

  void SetLowEnergyLimit(G4double e) {lowKinEnergy = e;};

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple*);
 
  G4bool IsInCharge(const G4ParticleDefinition*);

  virtual G4double ComputeDEDX(const G4Material*,
                       const G4ParticleDefinition*,
                             G4double kineticEnergy,
                             G4double cutEnergy);

  virtual G4double CrossSection(const G4Material*,
                        const G4ParticleDefinition*,
                              G4double kineticEnergy,
                              G4double cutEnergy,
                              G4double maxEnergy);

  virtual G4DynamicParticle* SampleSecondary(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);

  virtual G4std::vector<G4DynamicParticle*>* SampleSecondaries(
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

  G4double CrossSectionPerAtom(G4double Z,
                               G4double kineticEnergy,
                               G4double tmin,
                               G4double tmax,
                               G4double tmaxSecondary);

  G4double DifCrossSectionPerAtom(G4double kineticEnergy,
                                  G4double knockonEnergy,
                                  G4double tmaxSecondary);

  // hide assignment operator 
  G4MuBetheBlochModel & operator=(const  G4MuBetheBlochModel &right);
  G4MuBetheBlochModel(const  G4MuBetheBlochModel&);

  const G4ParticleDefinition* particle;
  G4double mass;
  G4double chargeSquare;
  G4double ratio;
  G4double highKinEnergy;
  G4double lowKinEnergy;
  G4double twoln10;
  G4double bg2lim;
  G4double taulim;
  G4double qc;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MuBetheBlochModel::MaxSecondaryEnergy(
          const G4ParticleDefinition* p,
                G4double kinEnergy) 
{

  G4double gamma= kinEnergy/mass + 1.0;
  G4double tmax = 2.0*electron_mass_c2*(gamma*gamma - 1.) /
                  (1. + 2.0*gamma*ratio + ratio*ratio);
  
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MuBetheBlochModel::MaxSecondaryEnergy(const G4DynamicParticle* dp)
{

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double gamma= kineticEnergy/mass + 1.0;
  G4double tmax = 2.0*electron_mass_c2*(gamma*gamma - 1.) /
                  (1. + 2.0*gamma*ratio + ratio*ratio);
  
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
