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
// $Id: G4eeToHadronsMultiModel.hh,v 1.1 2005/05/18 10:12:32 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eeToHadronsMultiModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 18.05.2005
//
// Modifications:
//

//
// Class Description: vector of e+e- -> hadrons models
//

// -------------------------------------------------------------------
//

#ifndef G4eeToHadronsMultiModel_h
#define G4eeToHadronsMultiModel_h 1

#include "G4VEmModel.hh"
#include "G4eeToHadronsModel.hh"
#include "Randomize.hh"
#include <vector>

class G4eeCrossSections;

class G4eeToHadronsMultiModel : public G4VEmModel
{

public:

  G4eeToHadronsMultiModel(G4int ver=0, const G4String& nam = "eeToHadrons");

  virtual ~G4eeToHadronsMultiModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double CrossSectionPerVolume(const G4Material*,
                                         const G4ParticleDefinition*,
                                         G4double kineticEnergy,
                                         G4double cutEnergy = 0.0,
                                         G4double maxEnergy = DBL_MAX);

  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin = 0.0,
                                      G4double maxEnergy = DBL_MAX);

  void PrintInfo();

  // Set the factor to artificially increase the crossSection (default 1)
  void SetCrossSecFactor(G4double fac);

private:

  // hide assignment operator
  G4eeToHadronsMultiModel & operator=(const  G4eeToHadronsMultiModel &right);
  G4eeToHadronsMultiModel(const  G4eeToHadronsMultiModel&);

  G4eeCrossSections*               cross;

  std::vector<G4eeToHadronsModel*> models;

  std::vector<G4double>            ekinMin;
  std::vector<G4double>            ekinPeak;
  std::vector<G4double>            ekinMax;
  std::vector<G4double>            cumSum;

  G4double                         thKineticEnergy;
  G4double                         maxKineticEnergy;
  G4double                         csFactor;

  G4int                            nModels;
  G4int                            verbose;
  G4bool                           isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToHadronsMultiModel::CrossSectionPerVolume(
                                 const G4Material* material,
                                 const G4ParticleDefinition*,
                                 G4double kineticEnergy,
                                 G4double, G4double)
{
  G4double res = 0.0;
  if (kineticEnergy > thKineticEnergy) {
    for(G4int i=0; i<nModels; i++) {
      if(kineticEnergy >= ekinMin[i] && kineticEnergy <= ekinMax[i])
        res += (models[i])->CrossSectionPerVolume(material,0,kineticEnergy);
      cumSum[i] = res;
    }
  }
  return res*csFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::vector<G4DynamicParticle*>* G4eeToHadronsMultiModel::SampleSecondaries(
                                const G4MaterialCutsCouple* couple,
                                const G4DynamicParticle* dp,
                                G4double, G4double)
{
  std::vector<G4DynamicParticle*>* newp = 0;
  G4double kinEnergy = dp->GetKineticEnergy();
  if (kinEnergy > thKineticEnergy) {
    G4double q = cumSum[nModels-1]*G4UniformRand();
    for(G4int i=0; i<nModels; i++) {
      if(q <= cumSum[i]) {
        newp = (models[i])->SampleSecondaries(couple,dp);
	break;
      }
    }
  }
  return newp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
