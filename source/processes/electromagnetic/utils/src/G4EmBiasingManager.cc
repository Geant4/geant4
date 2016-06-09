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
// $Id: G4EmBiasingManager.cc,v 1.88 2010-08-17 17:36:59 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
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
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmBiasingManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "Randomize.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmBiasingManager::G4EmBiasingManager() 
  : nForcedRegions(0),nSecBiasedRegions(0),
    currentStepLimit(0.0),startTracking(true)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmBiasingManager::~G4EmBiasingManager()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmBiasingManager::Initialise(const G4ParticleDefinition& part,
				    const G4String& procName, G4int verbose)
{
  //G4cout << "G4EmBiasingManager::Initialise for "
  //	 << part.GetParticleName()
  //	 << " and " << procName << G4endl;
  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  if(0 < nForcedRegions) { idxForcedCouple.resize(numOfCouples, -1); }
  if(0 < nSecBiasedRegions) { idxSecBiasedCouple.resize(numOfCouples, -1); }

  // Deexcitation
  for (size_t j=0; j<numOfCouples; ++j) {
    const G4MaterialCutsCouple* couple =
      theCoupleTable->GetMaterialCutsCouple(j);
    const G4ProductionCuts* pcuts = couple->GetProductionCuts();
    if(0 <  nForcedRegions) {
      for(G4int i=0; i<nForcedRegions; ++i) {
	if(forcedRegions[i]) {
	  if(pcuts == forcedRegions[i]->GetProductionCuts()) { 
	    idxForcedCouple[j] = i;
	    break; 
	  }
	}
      }
    }
    if(0 < nSecBiasedRegions) { 
      for(G4int i=0; i<nSecBiasedRegions; ++i) {
	if(secBiasedRegions[i]) {
	  if(pcuts == secBiasedRegions[i]->GetProductionCuts()) { 
	    idxSecBiasedCouple[j] = i;
	    break; 
	  }
	}
      }
    }
  }
  if (nForcedRegions > 0 && 0 < verbose) {
    G4cout << " Forced Interaction is activated for "
	   << part.GetParticleName() << " and " 
	   << procName 
	   << " inside G4Regions: " << G4endl;
    for (G4int i=0; i<nForcedRegions; ++i) {
      const G4Region* r = forcedRegions[i];
      if(r) { G4cout << "           " << r->GetName() << G4endl; }
    }
  }
  if (nSecBiasedRegions > 0 && 0 < verbose) {
    G4cout << " Secondary biasing is activated for " 
	   << part.GetParticleName() << " and " 
	   << procName 
	   << " inside G4Regions: " << G4endl;
    for (G4int i=0; i<nSecBiasedRegions; ++i) {
      const G4Region* r = secBiasedRegions[i];
      if(r) { 
	G4cout << "           " << r->GetName() 
	       << "  BiasingWeight= " << secBiasedWeight[i] << G4endl; 
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmBiasingManager::ActivateForcedInteraction(G4double val, 
						   const G4String& rname)
{
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  G4String name = rname;
  if(name == "" || name == "world" || name == "World") {
    name = "DefaultRegionForTheWorld";
  }
  const G4Region* reg = regionStore->GetRegion(name, false);
  if(!reg) { 
    G4cout << "### G4EmBiasingManager::ForcedInteraction WARNING: "
	   << " G4Region <"
	   << rname << "> is unknown" << G4endl;
    return; 
  }

  // the region is in the list
  if (0 < nForcedRegions) {
    for (G4int i=0; i<nForcedRegions; ++i) {
      if (reg == forcedRegions[i]) {
	lengthForRegion[i] = val; 
        return;
      }
    }
  }
  if(val < 0.0) { 
    G4cout << "### G4EmBiasingManager::ForcedInteraction WARNING: "
	   << val << " < 0.0, so no activation for the G4Region <"
	   << rname << ">" << G4endl;
    return; 
  }

  // new region 
  forcedRegions.push_back(reg);
  lengthForRegion.push_back(val);
  ++nForcedRegions;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4EmBiasingManager::ActivateSecondaryBiasing(const G4String& rname, 
					     G4double factor,
					     G4double energyLimit)
{
  //G4cout << "G4EmBiasingManager::ActivateSecondaryBiasing: "
  //	 << rname << " F= " << factor << " E(MeV)= " << energyLimit/MeV
  //	 << G4endl; 
  if(0.0 >= factor) { return; }
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  G4String name = rname;
  if(name == "" || name == "world" || name == "World") {
    name = "DefaultRegionForTheWorld";
  }
  const G4Region* reg = regionStore->GetRegion(name, false);
  if(!reg) { 
    G4cout << "### G4EmBiasingManager::ActivateBremsstrahlungSplitting WARNING: "
	   << " G4Region <"
	   << rname << "> is unknown" << G4endl;
    return; 
  }

  G4int nsplit = 0;
  G4double w = 1.0/factor;

  if(factor >= 1.0) {
    nsplit = G4int(factor + 0.5);
    w = 1.0/G4double(nsplit); 
  }

  // the region is in the list
  if (0 < nSecBiasedRegions) {
    for (G4int i=0; i<nSecBiasedRegions; ++i) {
      if (reg == secBiasedRegions[i]) {
	secBiasedWeight[i] = w;
        nBremSplitting[i]  = nsplit; 
	secBiasedEnegryLimit[i] = energyLimit;
        return;
      }
    }
  }
  if(1 == nsplit) { 
    G4cout << "### G4EmBiasingManager::ActivateSecondaryBiasing WARNING: "
	   << nsplit << " = 1, so no activation for the G4Region <"
	   << rname << ">" << G4endl;
    return; 
  }

  // new region 
  secBiasedRegions.push_back(reg);
  secBiasedWeight.push_back(w);
  nBremSplitting.push_back(nsplit);
  secBiasedEnegryLimit.push_back(energyLimit);
  ++nSecBiasedRegions;
  //G4cout << "nSecBiasedRegions= " << nSecBiasedRegions << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmBiasingManager::GetStepLimit(G4int coupleIdx, 
					  G4double previousStep)
{
  if(startTracking) {
    startTracking = false;
    G4int i = idxForcedCouple[coupleIdx];
    if(i < 0) {
      currentStepLimit = DBL_MAX;
    } else {
      currentStepLimit = lengthForRegion[i];
      if(currentStepLimit > 0.0) { currentStepLimit *= G4UniformRand(); }
    }
  } else {
    currentStepLimit -= previousStep;
  }
  if(currentStepLimit < 0.0) { currentStepLimit = 0.0; }
  return currentStepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EmBiasingManager::ApplySecondaryBiasing(std::vector<G4DynamicParticle*>& vd,
					  G4int coupleIdx)
{
  G4double weight = 1.0;
  size_t n = vd.size();
  G4int i = idxSecBiasedCouple[coupleIdx];
  if(0 <= i && 0 < n) {

    // apply biasing if first secondary has energy below the threshold
    if(vd[0]->GetKineticEnergy() < secBiasedEnegryLimit[i]) {
      weight = secBiasedWeight[i];
      G4int nsplit = nBremSplitting[i];

      // splitting
      if(1 < nsplit) {
	for(size_t k=0; k<n; ++k) {
	  const G4DynamicParticle* dp = vd[k];
	  for(G4int j=1; j<nsplit; ++j) {
	    G4DynamicParticle* dpnew = new G4DynamicParticle(*dp);
	    vd.push_back(dpnew);
	  }
	}
	// Russian roulette
      } else { 
	for(size_t k=0; k<n; ++k) {
	  const G4DynamicParticle* dp = vd[k];
	  if(G4UniformRand()*weight > 1.0) {
	    delete dp;
	    vd[k] = 0;
	  }
	}
      }
    }
  }
  return weight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4EmBiasingManager::ApplySecondaryBiasing(std::vector<G4Track*>& tr, 
					  G4double primaryWeight, 
					  G4int coupleIdx)
{
  G4double weight = primaryWeight;
  size_t n = tr.size();
  G4int i = idxSecBiasedCouple[coupleIdx];
  if(0 <= i && 0 < n) {

    weight *= secBiasedWeight[i];
    G4int nsplit = nBremSplitting[i];

    // splitting
    if(1 < nsplit) {
      for(size_t k=0; k<n; ++k) {
	G4Track* t = tr[k];
	t->SetWeight(weight);
	for(G4int j=1; j<nsplit; ++j) {
	  G4Track* tnew = new G4Track(*t);
          tnew->SetWeight(weight);
	  tr.push_back(tnew);
	}
      }
      // Russian roulette
    } else { 
      for(size_t k=0; k<n; ++k) {
	G4Track* t = tr[k];
        if(G4UniformRand()*secBiasedWeight[i] <= 1.0) {
          t->SetWeight(weight);
	} else {
	  delete t;
          tr[k] = 0;
	}
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
