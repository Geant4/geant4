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
// File name:     G4EmModelManager
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.05.2002
//
// Modifications:
//
// 03-12-02 V.Ivanchenko fix a bug in model selection
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 The set of models is defined for region (V.Ivanchenko)
//
// Class Description:
//
// It is the unified energy loss process it calculates the continuous
// energy loss for charged particles using a set of Energy Loss
// models valid for different energy regions. There are a possibility
// to create and access to dE/dx and range tables, or to calculate
// that information on fly.

// -------------------------------------------------------------------
//


#ifndef G4EmModelManager_h
#define G4EmModelManager_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4DataVector.hh"

class G4RegionModels
{

friend class G4EmModelManager;

public:

private:

  G4RegionModels(G4int nMod, G4std::vector<G4int>& list, G4DataVector& lowE);

  ~G4RegionModels();

  G4int SelectIndex(G4double e) const {
    G4int idx = 0;
    if (nModelsForRegion>1) {
      idx = nModelsForRegion;
      do {idx--;} while (idx && e <= lowKineticEnergy[idx]);
    }
    return theListOfModelIndexes[idx];
  };

  G4int ModelIndex(G4int n) const {
    return theListOfModelIndexes[n];
  };

  G4int NumberOfModels() const {
    return nModelsForRegion;
  };

  G4int              nModelsForRegion;
  G4int*             theListOfModelIndexes;
  G4double*          lowKineticEnergy;

};

#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4DynamicParticle.hh"

class G4Region;
class G4ParticleDefinition;
class G4DataVector;
class G4PhysicsVector;
class G4MaterialCutsCouple;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EmModelManager
{
public:

  G4EmModelManager();

 ~G4EmModelManager();

  void Clear();

  const G4DataVector* Initialise(const G4ParticleDefinition*,
                                 const G4ParticleDefinition*,
				       G4double,
                                       G4int);

  const G4DataVector* Cuts() const;

  const G4DataVector* SubCutoff() const;

  void FillDEDXVector(G4PhysicsVector*, const G4MaterialCutsCouple*);

  void FillLambdaVector(G4PhysicsVector*, const G4MaterialCutsCouple*);

  void FillSubLambdaVector(G4PhysicsVector*, const G4MaterialCutsCouple*);

  G4VEmModel* SelectModel(G4double& energy, size_t& index);

  G4double SampleFluctuations(const G4Material* material,
                              const G4DynamicParticle* dp,
                                    G4double& tmax,
		                    G4double& length,
		                    G4double& eloss,
				    G4double& preStepKinEnergy,
				    size_t&   index);

  void AddEmModel(G4int, G4VEmModel*, G4VEmFluctuationModel*, const G4Region*);

private:

  // hide  assignment operator

  G4EmModelManager(G4EmModelManager &);
  G4EmModelManager & operator=(const G4EmModelManager &right);

// =====================================================================

private:

  G4DataVector                theCuts;
  G4DataVector                theSubCuts;

  G4std::vector<G4VEmModel*>                models;
  G4std::vector<G4VEmFluctuationModel*>     flucModels;
  G4std::vector<const G4Region*>            regions;
  G4std::vector<G4int>                      orderOfModels;
  G4DataVector                              upperEkin;

  G4int                       nEmModels;
  G4int                       nRegions;
  G4int                       nCouples;

  G4std::vector<G4int>                      idxOfRegionModels;
  G4std::vector<G4RegionModels*>            setOfRegionModels;

  G4double                    minSubRange;

  const G4ParticleDefinition* particle;
  const G4ParticleDefinition* secondaryParticle;

  G4int                       verboseLevel;

  // cash
  G4int                       currentIdx;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4EmModelManager::SelectModel(G4double& kinEnergy, size_t& index)
{
  currentIdx = (setOfRegionModels[idxOfRegionModels[index]])->SelectIndex(kinEnergy);
  return models[currentIdx];
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4EmModelManager::SampleFluctuations(
                              const G4Material* material,
                              const G4DynamicParticle* dp,
                                    G4double& tmax,
		                    G4double& length,
		                    G4double& eloss,
				    G4double& kinEnergy,
				    size_t&   index)

{
  currentIdx = (setOfRegionModels[idxOfRegionModels[index]])->SelectIndex(kinEnergy);
  G4VEmFluctuationModel* fm = flucModels[currentIdx];
  G4double x = eloss;
  if(fm) x = fm->SampleFluctuations(material,dp,tmax,length,eloss);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4DataVector* G4EmModelManager::Cuts() const
{
  return &theCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4DataVector* G4EmModelManager::SubCutoff() const
{
  return &theSubCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

