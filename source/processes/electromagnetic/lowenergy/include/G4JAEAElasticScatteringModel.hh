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
/*
Authors:
M. Omer and R. Hajima  on   17 October 2016
contact:
omer.mohamed@jaea.go.jp and hajima.ryoichi@qst.go.jp
Publication Information:
1- M. Omer, R. Hajima, Including Delbrück scattering in Geant4,
Nucl. Instrum. Methods Phys. Res. Sect. B, vol. 405, 2017, pp. 43-49.,
https://doi.org/10.1016/j.nimb.2017.05.028
2- M. Omer, R. Hajima, Geant4 physics process for elastic scattering of gamma-rays,
JAEA Technical Report 2018-007, 2018.
https://doi.org/10.11484/jaea-data-code-2018-007
*/
//         on base of G4LivermoreRayleighModel
//

#ifndef G4JAEAElasticScatteringModel_h
#define G4JAEAElasticScatteringModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4ProductionCutsTable.hh"




class G4JAEAElasticScatteringModel : public G4VEmModel
{

public:

  G4JAEAElasticScatteringModel();

  virtual ~G4JAEAElasticScatteringModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual void InitialiseLocal(const G4ParticleDefinition*,
			       G4VEmModel* masterModel);

  virtual void InitialiseForElement(const G4ParticleDefinition*, G4int Z);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy,
                                      G4double Z,
                                      G4double A=0,
                                      G4double cut=0,
                                      G4double emax=DBL_MAX);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

  inline void SetLowEnergyThreshold(G4double);

    G4double distribution[181];
    G4double pdf[181];
    G4double cdf[181];

private:

  void ReadData(size_t Z, const char* path = 0);

  G4JAEAElasticScatteringModel & operator=(const G4JAEAElasticScatteringModel &right);
  G4JAEAElasticScatteringModel(const G4JAEAElasticScatteringModel&);

  G4bool isInitialised;
  G4int verboseLevel;

  G4double lowEnergyLimit;

  static G4int maxZ;
  static G4LPhysicsFreeVector* dataCS[101];

  G4ParticleChangeForGamma* fParticleChange;

};

inline void G4JAEAElasticScatteringModel::SetLowEnergyThreshold(G4double val)
{
  lowEnergyLimit = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
