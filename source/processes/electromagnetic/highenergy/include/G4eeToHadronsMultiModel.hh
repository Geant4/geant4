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
// $Id: G4eeToHadronsMultiModel.hh,v 1.7 2009-02-20 16:38:33 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4ParticleChangeForGamma.hh"
#include "G4TrackStatus.hh"
#include "Randomize.hh"
#include <vector>

class G4eeCrossSections;
class G4Vee2hadrons;

class G4eeToHadronsMultiModel : public G4VEmModel
{

public:

  G4eeToHadronsMultiModel(G4int ver=0, const G4String& nam = "eeToHadrons");

  virtual ~G4eeToHadronsMultiModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double CrossSectionPerVolume(const G4Material*,
					 const G4ParticleDefinition*,
					 G4double kineticEnergy,
					 G4double cutEnergy,
					 G4double maxEnergy);

  virtual G4double ComputeCrossSectionPerAtom(
                                         const G4ParticleDefinition*,
                                         G4double kineticEnergy,
                                         G4double Z, G4double A,
                                         G4double cutEnergy = 0.0,
                                         G4double maxEnergy = DBL_MAX);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin = 0.0,
				 G4double maxEnergy = DBL_MAX);

  virtual void PrintInfo();

  // Set the factor to artificially increase the crossSection (default 1)
  void SetCrossSecFactor(G4double fac);

  inline G4double ComputeCrossSectionPerElectron(
                                         const G4ParticleDefinition*,
                                         G4double kineticEnergy,
                                         G4double cutEnergy = 0.0,
                                         G4double maxEnergy = DBL_MAX);

private:

  void AddEEModel(G4Vee2hadrons*);

  // hide assignment operator
  G4eeToHadronsMultiModel & operator=(const  G4eeToHadronsMultiModel &right);
  G4eeToHadronsMultiModel(const  G4eeToHadronsMultiModel&);

  G4eeCrossSections*               cross;
  G4ParticleChangeForGamma*        fParticleChange;

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

inline G4double G4eeToHadronsMultiModel::ComputeCrossSectionPerElectron(
                                      const G4ParticleDefinition*,
				      G4double kineticEnergy,
				      G4double, G4double)
{
  G4double res = 0.0;
  if (kineticEnergy > thKineticEnergy) {
    for(G4int i=0; i<nModels; i++) {
      if(kineticEnergy >= ekinMin[i] && kineticEnergy <= ekinMax[i])
        res += (models[i])->ComputeCrossSectionPerElectron(0,kineticEnergy);
      cumSum[i] = res;
    }
  }
  return res*csFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
