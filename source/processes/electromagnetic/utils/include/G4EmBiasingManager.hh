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
// $Id: G4EmBiasingManager.hh,v 1.61 2010-08-17 17:36:59 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Region;
class G4DynamicParticle;
class G4Track;

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
  G4double ApplySecondaryBiasing(std::vector<G4DynamicParticle*>&, 
				 G4int coupleIdx);

  // Splitting or Russian roulette, G4Track may be deleted
  void ApplySecondaryBiasing(std::vector<G4Track*>&, G4double primaryWeight, 
			     G4int coupleIdx);

  inline G4bool SecondaryBiasingRegion(G4int coupleIdx);

  inline G4bool ForcedInteractionRegion(G4int coupleIdx);

  inline void ResetForcedInteraction();

private:

  // copy constructor and hide assignment operator
  G4EmBiasingManager(G4EmBiasingManager &);
  G4EmBiasingManager & operator=(const G4EmBiasingManager &right);

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

#endif
