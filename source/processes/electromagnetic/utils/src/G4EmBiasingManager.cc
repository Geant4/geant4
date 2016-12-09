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
// $Id: G4EmBiasingManager.cc 98860 2016-08-15 09:59:59Z gcosmo $
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
// 31-05-12 D. Sawkey put back in high energy limit for brem, russian roulette 
// 30-05-12 D. Sawkey  brem split gammas are unique; do weight tests for 
//          brem, russian roulette
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmBiasingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Track.hh"
#include "G4Electron.hh"
#include "G4VEmModel.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ParticleChangeForGamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmBiasingManager::G4EmBiasingManager() 
  : nForcedRegions(0),nSecBiasedRegions(0),eIonisation(nullptr),
    currentStepLimit(0.0),startTracking(true)
{
  fSafetyMin = 1.e-6*mm;
  theElectron = G4Electron::Electron();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmBiasingManager::~G4EmBiasingManager()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmBiasingManager::Initialise(const G4ParticleDefinition& part,
                                    const G4String& procName, G4int verbose)
{
  //G4cout << "G4EmBiasingManager::Initialise for "
  //         << part.GetParticleName()
  //         << " and " << procName << G4endl;
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
  //         << rname << " F= " << factor << " E(MeV)= " << energyLimit/MeV
  //         << G4endl; 
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  G4String name = rname;
  if(name == "" || name == "world" || name == "World") {
    name = "DefaultRegionForTheWorld";
  }
  const G4Region* reg = regionStore->GetRegion(name, false);
  if(!reg) { 
    G4cout << "### G4EmBiasingManager::ActivateBremsstrahlungSplitting "
           << "WARNING: G4Region <"
           << rname << "> is unknown" << G4endl;
    return; 
  }

  // Range cut
  G4int nsplit = 0;
  G4double w = factor;

  // splitting
  if(factor >= 1.0) {
    nsplit = G4lrint(factor);
    w = 1.0/G4double(nsplit);

    // Russian roulette 
  } else if(0.0 < factor) { 
    nsplit = 1;
    w = 1.0/factor;
  }

  // the region is in the list - overwrite parameters
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
  /*
    G4cout << "### G4EmBiasingManager::ActivateSecondaryBiasing: "
           << " nsplit= " << nsplit << " for the G4Region <"
           << rname << ">" << G4endl; 
  */

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
G4EmBiasingManager::ApplySecondaryBiasing(
                    std::vector<G4DynamicParticle*>& vd,
                    const G4Track& track,
                    G4VEmModel* currentModel,
                    G4ParticleChangeForLoss* pPartChange,
                    G4double& eloss,  
                    G4int coupleIdx,
                    G4double tcut, 
                    G4double safety)
{
  G4int index = idxSecBiasedCouple[coupleIdx];
  G4double weight = 1.0;
  if(0 <= index) {
    size_t n = vd.size();

    // the check cannot be applied per secondary particle
    // because weight correction is common, so the first
    // secondary is checked
    if(0 < n && vd[0]->GetKineticEnergy() < secBiasedEnegryLimit[index]) {

      G4int nsplit = nBremSplitting[index];

      // Range cut
      if(0 == nsplit) { 
        if(safety > fSafetyMin) { ApplyRangeCut(vd, track, eloss, safety); }

        // Russian Roulette
      } else if(1 == nsplit) { 
        weight = ApplyRussianRoulette(vd, index);

        // Splitting
      } else {
        G4double tmpEnergy = pPartChange->GetProposedKineticEnergy();
        G4ThreeVector tmpMomDir = pPartChange->GetProposedMomentumDirection();

        weight = ApplySplitting(vd, track, currentModel, index, tcut);

        pPartChange->SetProposedKineticEnergy(tmpEnergy);
        pPartChange->ProposeMomentumDirection(tmpMomDir);
      }
    }
  }
  return weight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EmBiasingManager::ApplySecondaryBiasing(
                  std::vector<G4DynamicParticle*>& vd,
                  const G4Track& track,
                  G4VEmModel* currentModel, 
                  G4ParticleChangeForGamma* pPartChange,
                  G4double& eloss,  
                  G4int coupleIdx,
                  G4double tcut, 
                  G4double safety)
{
  G4int index = idxSecBiasedCouple[coupleIdx];
  G4double weight = 1.0;
  if(0 <= index) {
    size_t n = vd.size();

    // the check cannot be applied per secondary particle
    // because weight correction is common, so the first
    // secondary is checked
    if(0 < n && vd[0]->GetKineticEnergy() < secBiasedEnegryLimit[index]) {

      G4int nsplit = nBremSplitting[index];

      // Range cut
      if(0 == nsplit) { 
        if(safety > fSafetyMin) { ApplyRangeCut(vd, track, eloss, safety); }

        // Russian Roulette
      } else if(1 == nsplit) { 
        weight = ApplyRussianRoulette(vd, index);

        // Splitting
      } else {
        G4double tmpEnergy = pPartChange->GetProposedKineticEnergy();
        G4ThreeVector tmpMomDir = pPartChange->GetProposedMomentumDirection();

        weight = ApplySplitting(vd, track, currentModel, index, tcut);

        pPartChange->SetProposedKineticEnergy(tmpEnergy);
        pPartChange->ProposeMomentumDirection(tmpMomDir);
      }
    }
  }
  return weight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4EmBiasingManager::ApplySecondaryBiasing(std::vector<G4Track*>& track,
                                          G4int coupleIdx)
{
  G4int index = idxSecBiasedCouple[coupleIdx];
  G4double weight = 1.0;
  if(0 <= index) {
    size_t n = track.size();

    // the check cannot be applied per secondary particle
    // because weight correction is common, so the first
    // secondary is checked
    if(0 < n && track[0]->GetKineticEnergy() < secBiasedEnegryLimit[index]) {

      G4int nsplit = nBremSplitting[index];

        // Russian Roulette only
      if(1 == nsplit) { 
        weight = secBiasedWeight[index];
        for(size_t k=0; k<n; ++k) {
          if(G4UniformRand()*weight > 1.0) {
            const G4Track* t = track[k];
            delete t;
            track[k] = 0;
          }
        }
      }
    }
  }
  return weight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4EmBiasingManager::ApplyRangeCut(std::vector<G4DynamicParticle*>& vd,
                                  const G4Track& track,
                                  G4double& eloss, G4double safety)
{
  size_t n = vd.size();
  if(!eIonisation) { 
    eIonisation = 
      G4LossTableManager::Instance()->GetEnergyLossProcess(theElectron);
  }
  if(eIonisation) { 
    for(size_t k=0; k<n; ++k) {
      const G4DynamicParticle* dp = vd[k];
      if(dp->GetDefinition() == theElectron) {
        G4double e = dp->GetKineticEnergy();
        if(eIonisation->GetRangeForLoss(e, track.GetMaterialCutsCouple()) 
           < safety) {
          eloss += e;
          delete dp;
          vd[k] = 0;
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4EmBiasingManager::ApplySplitting(std::vector<G4DynamicParticle*>& vd,
                                   const G4Track& track,
                                   G4VEmModel* currentModel, 
                                   G4int index,
                                   G4double tcut)
{
  // method is applied only if 1 secondary created PostStep 
  // in the case of many secodndaries there is a contrudition
  G4double weight = 1.0;
  size_t n = vd.size();
  G4double w = secBiasedWeight[index];

  if(1 != n || 1.0 <= w) { return weight; }

  G4double trackWeight = track.GetWeight();
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();

  G4int nsplit = nBremSplitting[index];

  // double splitting is supressed 
  if(1 < nsplit && trackWeight>w) {

    weight = w;
    // start from 1, because already one secondary created
    if(nsplit > (G4int)tmpSecondaries.size()) { 
      tmpSecondaries.reserve(nsplit);
    }
    const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
    for(G4int k=1; k<nsplit; ++k) {  
      tmpSecondaries.clear();
      currentModel->SampleSecondaries(&tmpSecondaries, couple, dynParticle, 
                                      tcut);
      for (size_t kk=0; kk<tmpSecondaries.size(); ++kk) {
        vd.push_back(tmpSecondaries[kk]);
      }
    }
  }
  return weight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
