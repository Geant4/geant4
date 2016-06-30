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
// $Id: G4EmBiasingManager.hh 95657 2016-02-17 13:03:36Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmBiasingManager
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 28.07.2011
//
// Modifications:
//
// Class Description:
//
// It is a class providing step limit for forced process biasing

// -------------------------------------------------------------------
//

#ifndef G4EmBiasingManager_h
#define G4EmBiasingManager_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Region;
class G4Track;
class G4VEnergyLossProcess;
class G4VEmModel;
class G4MaterialCutsCouple;
class G4ParticleChangeForLoss;
class G4ParticleChangeForGamma;

class G4EmBiasingManager 
{
public:

  G4EmBiasingManager();

  ~G4EmBiasingManager();

  void Initialise(const G4ParticleDefinition& part,
		  const G4String& procName, G4int verbose);

  // default parameters are possible
  void ActivateForcedInteraction(G4double length = 0.0, 
				 const G4String& r = "");

  // no default parameters
  void ActivateSecondaryBiasing(const G4String& region, G4double factor,
				G4double energyLimit);

  // return forced step limit
  G4double GetStepLimit(G4int coupleIdx, G4double previousStep);

  // return weight of splitting or Russian roulette
  // G4DynamicParticle may be deleted
  // two functions are required because of the different ParticleChange
  // ApplySecondaryBiasing() are wrappers

  // for G4VEmProcess 
  G4double ApplySecondaryBiasing(std::vector<G4DynamicParticle*>&, 
				 const G4Track& track,
				 G4VEmModel* currentModel,
				 G4ParticleChangeForGamma* pParticleChange,
				 G4double& eloss, 
   				 G4int coupleIdx,  
				 G4double tcut, 
				 G4double safety = 0.0);

  // for G4VEnergyLossProcess 
  G4double ApplySecondaryBiasing(std::vector<G4DynamicParticle*>&,
				 const G4Track& track,
				 G4VEmModel* currentModel,
				 G4ParticleChangeForLoss* pParticleChange,
				 G4double& eloss, 
   				 G4int coupleIdx,  
				 G4double tcut, 
				 G4double safety = 0.0);

  // for G4VEnergyLossProcess 
  G4double ApplySecondaryBiasing(std::vector<G4Track*>&,
				 G4int coupleIdx);

  inline G4bool SecondaryBiasingRegion(G4int coupleIdx);

  inline G4bool ForcedInteractionRegion(G4int coupleIdx);

  inline void ResetForcedInteraction();

private:

  void ApplyRangeCut(std::vector<G4DynamicParticle*>& vd,
		     const G4Track& track,
		     G4double& eloss, 
		     G4double safety);

  G4double ApplySplitting(std::vector<G4DynamicParticle*>& vd,
			  const G4Track& track,
			  G4VEmModel* currentModel, 
			  G4int index,
			  G4double tcut);

  inline G4double ApplyRussianRoulette(std::vector<G4DynamicParticle*>& vd,
				       G4int index); 

  // hide copy constructor and assignment operator
  G4EmBiasingManager(G4EmBiasingManager &) = delete;
  G4EmBiasingManager & operator=(const G4EmBiasingManager &right) = delete;

  G4int                        nForcedRegions;
  G4int                        nSecBiasedRegions;
  std::vector<const G4Region*> forcedRegions;
  std::vector<G4double>        lengthForRegion;
  std::vector<const G4Region*> secBiasedRegions;
  std::vector<G4double>        secBiasedWeight;
  std::vector<G4double>        secBiasedEnegryLimit;
  std::vector<G4int>           nBremSplitting;

  std::vector<G4int>           idxForcedCouple;
  std::vector<G4int>           idxSecBiasedCouple;

  std::vector<G4DynamicParticle*> tmpSecondaries;

  G4VEnergyLossProcess*         eIonisation;

  const G4ParticleDefinition*  theElectron;

  G4double fSafetyMin;
  G4double currentStepLimit;
  G4bool   startTracking;
};

inline G4bool 
G4EmBiasingManager::SecondaryBiasingRegion(G4int coupleIdx)
{
  G4bool res = false;
  if(nSecBiasedRegions > 0) {
    if(idxSecBiasedCouple[coupleIdx] >= 0) { res = true; }
  }
  return res;
}

inline G4bool G4EmBiasingManager::ForcedInteractionRegion(G4int coupleIdx)
{
  G4bool res = false;
  if(nForcedRegions > 0) {
    if(idxForcedCouple[coupleIdx] >= 0) { res = true; }
  }
  return res;
}

inline void G4EmBiasingManager::ResetForcedInteraction()
{
  startTracking = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double
G4EmBiasingManager::ApplyRussianRoulette(std::vector<G4DynamicParticle*>& vd,
					 G4int index)
{
  size_t n = vd.size();
  G4double weight = secBiasedWeight[index];
  for(size_t k=0; k<n; ++k) {
    if(G4UniformRand()*weight > 1.0) {
      const G4DynamicParticle* dp = vd[k];
      delete dp;
      vd[k] = nullptr;
    }
  }
  return weight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
