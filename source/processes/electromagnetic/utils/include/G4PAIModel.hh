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
// File name:     G4PAIModel
//
// Author:        V. Grichine based on Vladimir Ivanchenko  code
//
// Creation date: 05.10.2003
//
// Modifications:
//
//
// Class Description:
//
// Implementation of PAI model of energy loss and
// delta-electron production by heavy charged particles

// -------------------------------------------------------------------
//

#ifndef G4PAIModel_h
#define G4PAIModel_h 1

#include "G4VEmModel.hh"

class G4PAIModel : public G4VEmModel
{

public:

  G4PAIModel(const G4ParticleDefinition* p = 0, const G4String& nam = "PAI");

  ~G4PAIModel();

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

  std::vector<G4DynamicParticle*>* SampleSecondaries(
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
  G4PAIModel & operator=(const  G4PAIModel &right);
  G4PAIModel(const  G4PAIModel&);

  const G4ParticleDefinition* particle;
  G4double mass;
  G4double spin;
  G4double chargeSquare;
  G4double ratio;
  G4double highKinEnergy;
  G4double lowKinEnergy;
  G4double twoln10;
  G4double bg2lim; 
  G4double taulim;
  G4double qc;
};

/////////////////////////////////////////////////////////////////////

inline G4double G4PAIModel::MaxSecondaryEnergy(
          const G4ParticleDefinition*,
                G4double kinEnergy) 
{

  G4double gamma= kinEnergy/mass + 1.0;
  G4double tmax = 2.0*electron_mass_c2*(gamma*gamma - 1.) /
                  (1. + 2.0*gamma*ratio + ratio*ratio);
  
  return tmax;
}

/////////////////////////////////////////////////////////////////////////

inline G4double G4PAIModel::MaxSecondaryEnergy(const G4DynamicParticle* dp)
{

  G4double kineticEnergy = dp->GetKineticEnergy();
  G4double gamma= kineticEnergy/mass + 1.0;
  G4double tmax = 2.0*electron_mass_c2*(gamma*gamma - 1.) /
                  (1. + 2.0*gamma*ratio + ratio*ratio);
  
  return tmax;
}



#endif
