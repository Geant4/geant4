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
// Author: Sebastien Incerti
//         31 March 2012
//         on base of G4LivermoreRayleighModel
//

#ifndef G4LivermoreRayleighModel_h
#define G4LivermoreRayleighModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4ProductionCutsTable.hh"

class G4LivermoreRayleighModel : public G4VEmModel
{

public:

  G4LivermoreRayleighModel();

  virtual ~G4LivermoreRayleighModel();

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

private:

  void ReadData(size_t Z, const char* path = 0);

  G4LivermoreRayleighModel & operator=(const G4LivermoreRayleighModel &right);
  G4LivermoreRayleighModel(const G4LivermoreRayleighModel&);

  G4bool isInitialised;
  G4int verboseLevel;

  G4double lowEnergyLimit;  

  static G4int maxZ;
  static G4LPhysicsFreeVector* dataCS[101];

  G4ParticleChangeForGamma* fParticleChange;

};

inline void G4LivermoreRayleighModel::SetLowEnergyThreshold(G4double val)
{
  lowEnergyLimit = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
