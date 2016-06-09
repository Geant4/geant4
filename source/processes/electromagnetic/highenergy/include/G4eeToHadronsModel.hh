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
// $Id: G4eeToHadronsModel.hh,v 1.1 2004/11/19 18:44:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eeToHadronsModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 25.10.2003
//
// Modifications:
//

//
// Class Description:
//

// -------------------------------------------------------------------
//

#ifndef G4eeToHadronsModel_h
#define G4eeToHadronsModel_h 1

#include "G4VEmModel.hh"

class G4PhysicsVector;
class G4Vee2hadrons;

class G4eeToHadronsModel : public G4VEmModel
{

public:

  G4eeToHadronsModel(const G4Vee2hadrons*, G4int ver=0,
                     const G4String& nam = "eeToHadrons");

  virtual ~G4eeToHadronsModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double HighEnergyLimit(const G4ParticleDefinition* p = 0);

  G4double LowEnergyLimit(const G4ParticleDefinition* p = 0);

  G4double PeakEnergy() const;

  void SetHighEnergyLimit(G4double e);

  void SetLowEnergyLimit(G4double e);

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

  G4DynamicParticle* GenerateCMPhoton(G4double);

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                                    G4double kinEnergy);

private:

  void ComputeCMCrossSectionPerElectron();

  // hide assignment operator
  G4eeToHadronsModel & operator=(const  G4eeToHadronsModel &right);
  G4eeToHadronsModel(const  G4eeToHadronsModel&);

  const G4Vee2hadrons*  model;
  G4ParticleDefinition* theGamma;
  G4PhysicsVector*      crossPerElectron;
  G4PhysicsVector*      crossBornPerElectron;
  G4bool                isInitialised;
  G4int                 nbins;
  G4int                 verbose;

  G4double              lowKinEnergy;
  G4double              peakKinEnergy;
  G4double              highKinEnergy;

  G4double              emin;
  G4double              epeak;
  G4double              emax;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToHadronsModel::MaxSecondaryEnergy(
          const G4ParticleDefinition*,
                G4double kinEnergy)
{
  return kinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToHadronsModel::MaxSecondaryEnergy(const G4DynamicParticle* dp)
{
  return dp->GetKineticEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToHadronsModel::MinEnergyCut(const G4ParticleDefinition*,
                                                 const G4MaterialCutsCouple*)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToHadronsModel::PeakEnergy() const
{
  return peakKinEnergy;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
