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
// $Id: G4LowEnergyBraggModel.hh,v 1.6 2006/06/29 19:35:53 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4LowEnergyBraggModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 17.04.2004
//
// Modifications:
//
//
// Class Description:
//
// Implementation of energy loss and delta-electron production
// by heavy slow charged particles using eveluated data

// -------------------------------------------------------------------
//

#ifndef G4LowEnergyBraggModel_h
#define G4LowEnergyBraggModel_h 1

#include "G4VEmModel.hh"

class G4hParametrisedLossModel;

class G4LowEnergyBraggModel : public G4VEmModel
{

public:

  G4LowEnergyBraggModel(const G4ParticleDefinition* p = 0, const G4String& nam = "LEBragg");

  ~G4LowEnergyBraggModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double HighEnergyLimit(const G4ParticleDefinition*);

  G4double LowEnergyLimit(const G4ParticleDefinition*);

  void SetHighEnergyLimit(G4double e) {highKinEnergy = e;};

  void SetLowEnergyLimit(G4double e) {lowKinEnergy = e;};

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple*);

  G4bool IsInCharge(const G4ParticleDefinition*);

  G4double ComputeDEDX(const G4MaterialCutsCouple*,
                       const G4ParticleDefinition*,
                             G4double kineticEnergy,
                             G4double cutEnergy);

  G4double CrossSection(const G4MaterialCutsCouple*,
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

  void SetParameterization(G4hParametrisedLossModel* p);

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                                    G4double kinEnergy);

private:

  void SetParticle(const G4ParticleDefinition* p);

  // hide assignment operator
  G4LowEnergyBraggModel & operator=(const  G4LowEnergyBraggModel &right);
  G4LowEnergyBraggModel(const  G4LowEnergyBraggModel&);

  const G4ParticleDefinition* particle;
  G4hParametrisedLossModel* parameterization;
  G4double mass;
  G4double spin;
  G4double chargeSquare;
  G4double massRate;
  G4double ratio;
  G4double highKinEnergy;
  G4double lowKinEnergy;
  G4double lowestKinEnergy;

  G4bool   isIon;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4LowEnergyBraggModel::MaxSecondaryEnergy(
          const G4ParticleDefinition* pd,
                G4double kinEnergy)
{
  if(isIon) SetParticle(pd);
  G4double tau  = kinEnergy/mass;
  G4double tmax = 2.0*electron_mass_c2*tau*(tau + 2.) /
                  (1. + 2.0*(tau + 1.)*ratio + ratio*ratio);
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4LowEnergyBraggModel::MaxSecondaryEnergy(const G4DynamicParticle* dp)
{
  if(isIon) {
    mass =  dp->GetMass();
    ratio = electron_mass_c2/mass;
  }
  G4double tau  = dp->GetKineticEnergy()/mass;
  G4double tmax = 2.0*electron_mass_c2*tau*(tau + 2.) /
                  (1. + 2.0*(tau + 1.)*ratio + ratio*ratio);
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
