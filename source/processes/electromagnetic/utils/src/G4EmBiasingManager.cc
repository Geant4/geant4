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
#include "G4PhysicalConstants.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Track.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4VEmModel.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmBiasingManager::G4EmBiasingManager()
  : fDirectionalSplittingTarget(0.0,0.0,0.0)
{
  fSafetyMin = 1.e-6*mm;
  theElectron = G4Electron::Electron();
  theGamma    = G4Gamma::Gamma();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmBiasingManager::~G4EmBiasingManager() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmBiasingManager::Initialise(const G4ParticleDefinition& part,
                                    const G4String& procName, G4int verbose)
{
  //G4cout << "G4EmBiasingManager::Initialise for "
  //         << part.GetParticleName()
  //         << " and " << procName << G4endl;
  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = (G4int)theCoupleTable->GetTableSize();

  if(0 < nForcedRegions) { idxForcedCouple.resize(numOfCouples, -1); }
  if(0 < nSecBiasedRegions) { idxSecBiasedCouple.resize(numOfCouples, -1); }

  // Deexcitation
  for (G4int j=0; j<numOfCouples; ++j) {
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

  G4EmParameters* param = G4EmParameters::Instance();
  SetDirectionalSplitting(param->GetDirectionalSplitting());
  if (fDirectionalSplitting) {
    SetDirectionalSplittingTarget(param->GetDirectionalSplittingTarget());
    SetDirectionalSplittingRadius(param->GetDirectionalSplittingRadius());
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
    if (fDirectionalSplitting) {
      G4cout << "     Directional splitting activated, with target position: "
             << fDirectionalSplittingTarget/cm
             << " cm; radius: "
             << fDirectionalSplittingRadius/cm
             << "cm." << G4endl;
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
  G4double weight = 1.;
  if(0 <= index) {
    std::size_t n = vd.size();

    // the check cannot be applied per secondary particle
    // because weight correction is common, so the first
    // secondary is checked
    if((0 < n && vd[0]->GetKineticEnergy() < secBiasedEnegryLimit[index])
          || fDirectionalSplitting) {

      G4int nsplit = nBremSplitting[index];

      // Range cut
      if(0 == nsplit) { 
        if(safety > fSafetyMin) { ApplyRangeCut(vd, track, eloss, safety); }

        // Russian Roulette
      } else if(1 == nsplit) { 
        weight = ApplyRussianRoulette(vd, index);

        // Splitting
      } else {
        if (fDirectionalSplitting) {
          weight = ApplyDirectionalSplitting(vd, track, currentModel, index, tcut);
        } else {
          G4double tmpEnergy = pPartChange->GetProposedKineticEnergy();
          G4ThreeVector tmpMomDir = pPartChange->GetProposedMomentumDirection();

          weight = ApplySplitting(vd, track, currentModel, index, tcut);

          pPartChange->SetProposedKineticEnergy(tmpEnergy);
          pPartChange->ProposeMomentumDirection(tmpMomDir);
        }
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
  G4double weight = 1.;
  if(0 <= index) {
    std::size_t n = vd.size();

    // the check cannot be applied per secondary particle
    // because weight correction is common, so the first
    // secondary is checked
    if((0 < n && vd[0]->GetKineticEnergy() < secBiasedEnegryLimit[index])
          || fDirectionalSplitting) {

      G4int nsplit = nBremSplitting[index];

      // Range cut
      if(0 == nsplit) { 
        if(safety > fSafetyMin) { ApplyRangeCut(vd, track, eloss, safety); }

        // Russian Roulette
      } else if(1 == nsplit) { 
        weight = ApplyRussianRoulette(vd, index);

        // Splitting
      } else {
        if (fDirectionalSplitting) {
          weight = ApplyDirectionalSplitting(vd, track, currentModel,
                                    index, tcut, pPartChange);
        } else {
          G4double tmpEnergy = pPartChange->GetProposedKineticEnergy();
          G4ThreeVector tmpMomDir = pPartChange->GetProposedMomentumDirection();

          weight = ApplySplitting(vd, track, currentModel, index, tcut);

          pPartChange->SetProposedKineticEnergy(tmpEnergy);
          pPartChange->ProposeMomentumDirection(tmpMomDir);
        }
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
  G4double weight = 1.;
  if(0 <= index) {
    std::size_t n = track.size();

    // the check cannot be applied per secondary particle
    // because weight correction is common, so the first
    // secondary is checked
    if(0 < n && track[0]->GetKineticEnergy() < secBiasedEnegryLimit[index]) {

      G4int nsplit = nBremSplitting[index];

        // Russian Roulette only
      if(1 == nsplit) { 
        weight = secBiasedWeight[index];
        for(std::size_t k=0; k<n; ++k) {
          if(G4UniformRand()*weight > 1.0) {
            const G4Track* t = track[k];
            delete t;
            track[k] = nullptr;
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
  std::size_t n = vd.size();
  if(!eIonisation) { 
    eIonisation = 
      G4LossTableManager::Instance()->GetEnergyLossProcess(theElectron);
  }
  if(eIonisation) { 
    for(std::size_t k=0; k<n; ++k) {
      const G4DynamicParticle* dp = vd[k];
      if(dp->GetDefinition() == theElectron) {
        G4double e = dp->GetKineticEnergy();
        if(eIonisation->GetRange(e, track.GetMaterialCutsCouple()) < safety) {
          eloss += e;
          delete dp;
          vd[k] = nullptr;
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4EmBiasingManager::CheckDirection(G4ThreeVector pos,
                                          G4ThreeVector momdir) const
{
  G4ThreeVector delta = fDirectionalSplittingTarget - pos;
  G4double angle = momdir.angle(delta);
  G4double dist = delta.cross(momdir).mag();
  if (dist <= fDirectionalSplittingRadius && angle < halfpi) {
    return true;
  }
  return false;
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
  // in the case of many secondaries there is a contradiction
  G4double weight = 1.;
  std::size_t n = vd.size();
  G4double w = secBiasedWeight[index];

  if(1 != n || 1.0 <= w) { return weight; }

  G4double trackWeight = track.GetWeight();
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();

  G4int nsplit = nBremSplitting[index];

  // double splitting is suppressed 
  if(1 < nsplit && trackWeight>w) {

    weight = w;
    if(nsplit > (G4int)tmpSecondaries.size()) { 
      tmpSecondaries.reserve(nsplit);
    }
    const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
    // start from 1, because already one secondary created
    for(G4int k=1; k<nsplit; ++k) {
      tmpSecondaries.clear();
      currentModel->SampleSecondaries(&tmpSecondaries, couple, dynParticle, 
                                      tcut);
      for (std::size_t kk=0; kk<tmpSecondaries.size(); ++kk) {
        vd.push_back(tmpSecondaries[kk]);
      }
    }
  }
  return weight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4EmBiasingManager::ApplyDirectionalSplitting(
                                   std::vector<G4DynamicParticle*>& vd,
                                   const G4Track& track,
                                   G4VEmModel* currentModel,
                                   G4int index,
                                   G4double tcut,
                                   G4ParticleChangeForGamma* partChange)
{
  // primary is gamma. do splitting/RR as appropriate
  // method applied for any number of secondaries

  G4double weight = 1.0;
  G4double w = secBiasedWeight[index];

  fDirectionalSplittingWeights.clear();
  if(1.0 <= w) {
    fDirectionalSplittingWeights.push_back(weight);
    return weight;
  }

  G4double trackWeight = track.GetWeight();
  G4int nsplit = nBremSplitting[index];

  // double splitting is suppressed
  if(1 < nsplit && trackWeight>w) {

    weight = w;
    const G4ThreeVector pos = track.GetPosition();

    G4bool foundPrimaryParticle = false;
    G4double primaryEnergy = 0.;
    G4ThreeVector primaryMomdir(0.,0.,0.);
    G4double primaryWeight = trackWeight;

    tmpSecondaries = vd;
    vd.clear();
    vd.reserve(nsplit);
    for (G4int k=0; k<nsplit; ++k) {
      if (k>0) {  // for k==0, SampleSecondaries has already been called
        tmpSecondaries.clear();
        // SampleSecondaries modifies primary info stored in partChange
        currentModel->SampleSecondaries(&tmpSecondaries,
                                        track.GetMaterialCutsCouple(),
                                        track.GetDynamicParticle(), tcut);
      }
      for (std::size_t kk=0; kk<tmpSecondaries.size(); ++kk) {
        if (tmpSecondaries[kk]->GetParticleDefinition() == theGamma) {
          if (CheckDirection(pos, tmpSecondaries[kk]->GetMomentumDirection())){
            vd.push_back(tmpSecondaries[kk]);
            fDirectionalSplittingWeights.push_back(1.);
          } else if (G4UniformRand() < w) {
            vd.push_back(tmpSecondaries[kk]);
            fDirectionalSplittingWeights.push_back(1./weight);
          } else {
            delete tmpSecondaries[kk];
            tmpSecondaries[kk] = nullptr;
          }
        } else if (k==0) { // keep charged 2ry from first splitting
          vd.push_back(tmpSecondaries[kk]);
          fDirectionalSplittingWeights.push_back(1./weight);
        } else {
          delete tmpSecondaries[kk];
          tmpSecondaries[kk] = nullptr;
        }
      }

      // primary
      G4double en = partChange->GetProposedKineticEnergy();
      if (en>0.) { // don't add if kinetic energy = 0
        G4ThreeVector momdir = partChange->GetProposedMomentumDirection();
        if (CheckDirection(pos,momdir)) {
          // keep only one primary; others are secondaries
          if (!foundPrimaryParticle) {
            primaryEnergy = en;
            primaryMomdir = momdir;
            foundPrimaryParticle = true;
            primaryWeight = weight;
          } else {
            auto dp = new G4DynamicParticle(theGamma,
                          partChange->GetProposedMomentumDirection(),
                          partChange->GetProposedKineticEnergy());
            vd.push_back(dp);
            fDirectionalSplittingWeights.push_back(1.);
          }
        } else if (G4UniformRand()<w) { // not going to target. play RR.
          if (!foundPrimaryParticle) {
            foundPrimaryParticle = true;
            primaryEnergy = en;
            primaryMomdir = momdir;
            primaryWeight = 1.;
          } else {
            auto dp = new G4DynamicParticle(theGamma,
                          partChange->GetProposedMomentumDirection(),
                          partChange->GetProposedKineticEnergy());
            vd.push_back(dp);
            fDirectionalSplittingWeights.push_back(1./weight);
          }
        }
      }
    }  // end of loop over nsplit

    partChange->ProposeWeight(primaryWeight);
    partChange->SetProposedKineticEnergy(primaryEnergy);
    partChange->ProposeMomentumDirection(primaryMomdir);
  } else {
    for (std::size_t i = 0; i < vd.size(); ++i) {
      fDirectionalSplittingWeights.push_back(1.);
    }
  }

  return weight;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmBiasingManager::GetWeight(G4int i)
{
  // normally return 1. If a directionally split particle survives RR,
  //  return 1./(splitting factor)
  if (fDirectionalSplittingWeights.size() >= (unsigned int)(i+1) ) {
    G4double w = fDirectionalSplittingWeights[i];
    fDirectionalSplittingWeights[i] = 1.; // ensure it's not used again
    return w;
  } else {
    return 1.;
  }
}

G4double
G4EmBiasingManager::ApplyDirectionalSplitting(
                                  std::vector<G4DynamicParticle*>& vd,
                                  const G4Track& track,
                                  G4VEmModel* currentModel,
                                  G4int index,
                                  G4double tcut)
{
  // primary is not a gamma. Do nothing with primary

  G4double weight = 1.0;
  G4double w = secBiasedWeight[index];

  fDirectionalSplittingWeights.clear();
  if(1.0 <= w) {
    fDirectionalSplittingWeights.push_back(weight);
    return weight;
  }

  G4double trackWeight = track.GetWeight();
  G4int nsplit = nBremSplitting[index];

  // double splitting is suppressed
  if(1 < nsplit && trackWeight>w) {

    weight = w;
    const G4ThreeVector pos = track.GetPosition();

    tmpSecondaries = vd;
    vd.clear();
    vd.reserve(nsplit);
    for (G4int k=0; k<nsplit; ++k) {
      if (k>0) {
        tmpSecondaries.clear();
        currentModel->SampleSecondaries(&tmpSecondaries,
                                        track.GetMaterialCutsCouple(),
                                        track.GetDynamicParticle(), tcut);
      }
      //for (auto sec : tmpSecondaries) {
      for (std::size_t kk=0; kk < tmpSecondaries.size(); ++kk) {
        if (CheckDirection(pos, tmpSecondaries[kk]->GetMomentumDirection())) {
          vd.push_back(tmpSecondaries[kk]);
          fDirectionalSplittingWeights.push_back(1.);
        } else if (G4UniformRand()<w) {
          vd.push_back(tmpSecondaries[kk]);
          fDirectionalSplittingWeights.push_back(1./weight);
        } else {
          delete tmpSecondaries[kk];
          tmpSecondaries[kk] = nullptr;
        }
      }
    }  // end of loop over nsplit
  } else { // no splitting was done; still need weights
    for (std::size_t i = 0; i < vd.size(); ++i) {
      fDirectionalSplittingWeights.push_back(1.0);
    }
  }
  return weight;
}
