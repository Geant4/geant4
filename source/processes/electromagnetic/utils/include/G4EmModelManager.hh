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
// Modifications: 03.12.2002 V.Ivanchenko fix a bug in model selection
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

#include "globals.hh"
#include "G4Material.hh"
#include "G4Track.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"

class G4Step;
class G4ParticleDefinition;
class G4VEmModel;
class G4DataVector;
class G4VParticleChange;

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

  void FillDEDXVector(G4PhysicsVector*, const G4Material*);

  void FillLambdaVector(G4PhysicsVector*, const G4Material*);

  void FillSubLambdaVector(G4PhysicsVector*, const G4Material*);

  G4VEmModel* SelectModel(G4double energy);

  void AddEmModel(G4VEmModel*, G4int);

private:

  // hide  assignment operator 

  G4EmModelManager(G4EmModelManager &);
  G4EmModelManager & operator=(const G4EmModelManager &right);

// =====================================================================

private:

  G4DataVector                theCuts;
  G4DataVector                theSubCuts;

  G4VEmModel*                 emModels[5];
  G4int                       nEmModels;
  G4int                       nmax;

  G4int                       order[5];
  G4double                    upperEkin[5];
  G4bool                      orderIsChanged;

  G4double                    minSubRange;

  const G4ParticleDefinition* particle;
  const G4ParticleDefinition* secondaryParticle;

  G4int                       verboseLevel;  

  // cash
  G4VEmModel*                 currentModel;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VEmModel* G4EmModelManager::SelectModel(G4double energy)
{
  if(nEmModels == 4 && energy > upperEkin[2]) {
    currentModel = emModels[3];
  } else if(nEmModels >= 3 && energy > upperEkin[1]) {
    currentModel = emModels[2];
  } else if(nEmModels >= 2 && energy > upperEkin[0]) {
    currentModel = emModels[1];
  } else {    
    currentModel = emModels[0];
  }  
  return currentModel;
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
 


