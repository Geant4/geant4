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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
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
// 26-03-03 Add GetDEDXDispersion (V.Ivanchenko)
// 13-04-03 Add startFromNull (V.Ivanchenko)
// 13-05-03 Add calculation of precise range (V.Ivanchenko)
// 21-07-03 Add UpdateEmModel method (V.Ivanchenko)
// 03-11-03 Substitute STL vector for G4RegionModels (V.Ivanchenko)
// 11-04-05 Remove access to fluctuation models (V.Ivanchenko)
// 10-01-06 PreciseRange -> CSDARange (V.Ivantchenko)
// 20-01-06 Introduce G4EmTableType and reducing number of methods (VI)
// 13-05-06 Add GetModel by index method (VI)
// 15-03-07 Add maxCutInRange (V.Ivanchenko)
// 08-04-08 Simplify Select method for only one G4RegionModel (VI)
// 03-08-09 Removed unused members and simplify model search if only one
//          model is used (VI)
// 14-07-11 Use pointer to the vector of cuts and not local copy (VI)
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
#include "G4EmTableType.hh"
#include "G4EmProcessSubType.hh"
#include "G4Region.hh"

#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4DynamicParticle.hh"
#include <iostream>

class G4RegionModels
{

friend class G4EmModelManager;

private:

  G4RegionModels(G4int nMod, std::vector<G4int>& indx, 
                 G4DataVector& lowE, const G4Region* reg);

  ~G4RegionModels();

  inline G4int SelectIndex(G4double e) const {
    G4int idx = 0;
    if (nModelsForRegion>1) {
      idx = nModelsForRegion;
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
      do {--idx;} while (idx > 0 && e <= lowKineticEnergy[idx]);
    }
    return theListOfModelIndexes[idx];
  };

  inline G4int ModelIndex(G4int n) const {
    return theListOfModelIndexes[n];
  };

  inline G4int NumberOfModels() const {
    return nModelsForRegion;
  };

  inline G4double LowEdgeEnergy(G4int n) const {
    return lowKineticEnergy[n];
  };

  inline const G4Region* Region() const {
    return theRegion;
  };

  G4RegionModels(G4RegionModels &) = delete;
  G4RegionModels & operator=(const G4RegionModels &right) = delete;

  const G4Region*    theRegion;
  G4int              nModelsForRegion;
  G4int*             theListOfModelIndexes;
  G4double*          lowKineticEnergy;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Region;
class G4ParticleDefinition;
class G4PhysicsVector;
class G4MaterialCutsCouple;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EmModelManager
{
public:

  G4EmModelManager();

  ~G4EmModelManager();

  void Clear();

  const G4DataVector* Initialise(const G4ParticleDefinition* part,
                                 const G4ParticleDefinition* secPart,
                                 G4int verb);

  void FillDEDXVector(G4PhysicsVector*, const G4MaterialCutsCouple*, 
                      G4EmTableType t = fRestricted);

  void FillLambdaVector(G4PhysicsVector*, const G4MaterialCutsCouple*, 
                        G4bool startFromNull = true, 
                        G4EmTableType t = fRestricted);

  void AddEmModel(G4int, G4VEmModel*, G4VEmFluctuationModel* fm, 
                  const G4Region* r);

  // Get model pointer from the model list
  G4VEmModel* GetModel(G4int idx, G4bool ver = false) const;

  // Get model pointer from the model list for a given material cuts couple
  // no check on material cuts couple index
  G4VEmModel* GetRegionModel(G4int idx, std::size_t index_couple);

  // total number of models for material cut couples
  // no check on material cuts couple index
  G4int NumberOfRegionModels(std::size_t index_couple) const;

  // Automatic documentation
  void DumpModelList(std::ostream& out, G4int verb);

  // Select model for given material cuts couple index
  inline G4VEmModel* SelectModel(G4double energy, std::size_t index);

  // Access to cuts
  inline const G4DataVector* Cuts() const;

  // Set flag of fluorescence
  inline void SetFluoFlag(G4bool val);

  // total number of models
  inline G4int NumberOfModels() const;

  // hide  assignment operator
  G4EmModelManager(G4EmModelManager &) = delete;
  G4EmModelManager & operator=(const G4EmModelManager &right) = delete;

private:

  const G4ParticleDefinition* particle = nullptr;
  const G4DataVector*         theCuts = nullptr;
  G4DataVector*               theCutsNew = nullptr;

  // may be changed in run time
  G4RegionModels*             currRegionModel = nullptr;
  G4VEmModel*                 currModel = nullptr;

  G4int                       nEmModels = 0;
  G4int                       nRegions = 0;

  G4int                       verboseLevel = 0;
  G4bool                      severalModels = true;
  G4bool                      fluoFlag = false;

  std::vector<G4VEmModel*>             models;
  std::vector<G4VEmFluctuationModel*>  flucModels;
  std::vector<const G4Region*>         regions;
  std::vector<G4int>                   orderOfModels;
  std::vector<G4int>                   isUsed;

  std::vector<G4int>            idxOfRegionModels;
  std::vector<G4RegionModels*>  setOfRegionModels;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4VEmModel* G4EmModelManager::SelectModel(G4double kinEnergy, std::size_t index)
{
  if(severalModels) {
    if(nRegions > 1) {
      currRegionModel = setOfRegionModels[idxOfRegionModels[index]];
    }
    currModel = models[currRegionModel->SelectIndex(kinEnergy)];
  }
  return currModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4DataVector* G4EmModelManager::Cuts() const
{
  return theCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4EmModelManager::SetFluoFlag(G4bool val)
{
  fluoFlag = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4EmModelManager::NumberOfModels() const
{
  return nEmModels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

