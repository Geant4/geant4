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
// $Id: G4EmCalculator.cc,v 1.3 2004-07-22 14:44:07 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4EmCalculator
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 28.06.2004
//
// Modifications:
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmCalculator.hh"
#include "G4LossTableManager.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VMultipleScattering.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ionEffectiveCharge.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmCalculator::G4EmCalculator()
{
  manager = G4LossTableManager::Instance();
  nLocalMaterials    = 0;
  currentCoupleIndex = 0;
  currentCouple      = 0;
  currentMaterial    = 0;
  currentParticle    = 0;
  currentLambda      = 0;
  chargeSquare       = 1.0;
  massRatio          = 1.0;
  ionEffCharge       = new G4ionEffectiveCharge();
  isIon              = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmCalculator::~G4EmCalculator()
{
  delete ionEffCharge;
  for (G4int i=0; i<nLocalMaterials; i++) {
    delete localCouples[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetDEDX(const G4ParticleDefinition* p, const G4Material* mat,
                                       G4double kinEnergy)
{
  G4double res = 0.0;
  FindMaterial(mat);
  UpdateParticle(p, kinEnergy);
  if(currentCouple) res = manager->GetDEDX(p, kinEnergy, currentCouple);
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetRange(const G4ParticleDefinition* p, const G4Material* mat,
                                        G4double kinEnergy)
{
  G4double res = 0.0;
  FindMaterial(mat);
  UpdateParticle(p, kinEnergy);
  if(currentCouple) res = manager->GetRange(p, kinEnergy, currentCouple);
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetKinEnergy(const G4ParticleDefinition* p, 
                                      const G4Material* mat,
                                            G4double range)
{
  G4double res = 0.0;
  FindMaterial(mat);
  UpdateParticle(p, 1.0*GeV);
  if(currentCouple) res = manager->GetEnergy(p, range, currentCouple);
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetCrossSection(const G4ParticleDefinition* p, 
                                         const G4Material* mat,
                                         const G4String& processName, G4double kinEnergy)
{
  G4double res = 0.0;
  FindMaterial(mat);

  if(currentCouple) {
    FindLambdaTable(p, processName);
    UpdateParticle(p, kinEnergy);
    if(currentLambda) {
      G4bool b;
      G4double e = kinEnergy*massRatio;
      res = (((*currentLambda)[currentCoupleIndex])->GetValue(e,b))*chargeSquare;
      res /= currentMaterial->GetTotNbOfAtomsPerVolume();
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetMinFreePath(const G4ParticleDefinition* p, 
                                        const G4Material* mat,
                                        const G4String& processName, 
                                              G4double kinEnergy)
{
  G4double res = DBL_MAX;
  G4double x = GetCrossSection(p, mat, processName, kinEnergy);
  if(x > 0.0) res = 1.0/(x*(currentMaterial->GetTotNbOfAtomsPerVolume()));
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::PrintDEDXTable(const G4ParticleDefinition* p)
{
  const G4VEnergyLossProcess* elp = FindEnergyLossProcess(p);
  G4cout << "##### DEDX Table for " << p->GetParticleName() << G4endl;
  if(elp) G4cout << *(elp->DEDXTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::PrintRangeTable(const G4ParticleDefinition* p)
{
  const G4VEnergyLossProcess* elp = FindEnergyLossProcess(p);
  G4cout << "##### Range Table for " << p->GetParticleName() << G4endl;
  if(elp) G4cout << *(elp->RangeTableForLoss()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::PrintInverseRangeTable(const G4ParticleDefinition* p)
{
  const G4VEnergyLossProcess* elp = FindEnergyLossProcess(p);
  G4cout << "##### Inverse Range Table for " << p->GetParticleName() << G4endl;
  if(elp) G4cout << *(elp->InverseRangeTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeDEDX(const G4ParticleDefinition* p, 
                                     const G4Material* mat,
                                     const G4String& processName,
                                           G4double kinEnergy, 
                                           G4double cut, 
                                           size_t idxRegion)
{
  G4double res = 0.0;
  FindCouple(mat, cut);
  UpdateParticle(p, kinEnergy);
  FindEmModel(p, processName, kinEnergy, idxRegion);
  if(currentModel) {
    G4double e = kinEnergy*massRatio;
    res = currentModel->ComputeDEDX(currentCouple, baseParticle, e, cut)*chargeSquare;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeCrossSection(const G4ParticleDefinition* p,
                                             const G4Material* mat,
                                             const G4String& processName,
                                                   G4double kinEnergy, 
                                                   G4double cut, 
                                                   size_t idxRegion)
{
  G4double res = 0.0;
  FindCouple(mat, cut);
  UpdateParticle(p, kinEnergy);
  FindEmModel(p, processName, kinEnergy, idxRegion);
  if(currentModel) {
    G4double e = kinEnergy*massRatio;
    res = currentModel->CrossSection(currentCouple, baseParticle, e, cut, e)*chargeSquare;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeMinFreePath(const G4ParticleDefinition* p, 
                                            const G4Material* mat,
                                            const G4String& processName,
                                                  G4double kinEnergy, 
                                                  G4double cut, 
                                                  size_t idxRegion)
{
  G4double mfp = DBL_MAX;
  G4double x = ComputeCrossSection(p, mat, processName, kinEnergy, cut, idxRegion);
  if(x > 0.0) mfp = 1.0/(x*(currentMaterial->GetTotNbOfAtomsPerVolume()));
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::UpdateParticle(const G4ParticleDefinition* p, G4double kinEnergy)
{
  if(p != currentParticle) {
    currentParticle = p;
    const G4VEnergyLossProcess* proc = FindEnergyLossProcess(p);
    baseParticle = proc->BaseParticle();
    if(p->GetParticleType() == "nucleus") {
      isIon = true;
    } else {
      isIon = false;
      G4double q = p->GetPDGCharge()/eplus;
      chargeSquare = q*q;  
    }
    massRatio = 1.0;
    if(baseParticle) massRatio = baseParticle->GetPDGMass()/p->GetPDGMass();
  }
  if(isIon)
    chargeSquare = ionEffCharge->EffectiveChargeSquareRatio(p, currentMaterial, kinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::FindMaterial(const G4Material* material)
{
  if(currentMaterial == material) return;

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  // Search for the first couple with the given material
  currentCouple = 0;
  const G4MaterialCutsCouple* cpl    = 0;
  for(size_t i=0; i<numOfCouples; i++) {
    cpl = theCoupleTable->GetMaterialCutsCouple(i);
    currentMaterial = cpl->GetMaterial();
    if(currentMaterial == material) {
      currentCouple = cpl;
      currentCoupleIndex = i;
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::FindCouple(const G4Material* material, G4double cut)
{
  currentMaterial = material;
  for (G4int i=0; i<nLocalMaterials; i++) {
    if(material == localMaterials[i] && cut == localCuts[i]) {
      currentCouple = localCouples[i];
      currentCoupleIndex = currentCouple->GetIndex();
      currentCut = cut;
      return;
    }
  }
  const G4MaterialCutsCouple* cc = new G4MaterialCutsCouple(material);
  localMaterials.push_back(material);
  localCouples.push_back(cc);
  localCuts.push_back(cut);
  nLocalMaterials++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::FindLambdaTable(const G4ParticleDefinition* p,
                                     const G4String& processName)
{
  // Search for the process
  if (p != currentParticle || processName != currentName) {
    currentName     = processName;
    currentLambda   = 0;
    G4LossTableManager* manager = G4LossTableManager::Instance();
    const std::vector<G4VEnergyLossProcess*> vel = manager->GetEnergyLossProcessVector();
    G4int n = vel.size();
    for(G4int i=0; i<n; i++) {
      if((vel[i])->GetProcessName() == currentName && (vel[i])->Particle() == p) {
        currentLambda = (vel[i])->LambdaTable();
	break;
      }
    }
    if(!currentLambda) {
      const std::vector<G4VEmProcess*> vem = manager->GetEmProcessVector();
      G4int n = vem.size();
      for(G4int i=0; i<n; i++) {
        if((vem[i])->GetProcessName() == currentName && (vem[i])->Particle() == p) {
          currentLambda = (vem[i])->LambdaTable();
	  break;
        }
      }
    }
    if(!currentLambda) {
      const std::vector<G4VMultipleScattering*> vmsc = manager->GetMultipleScatteringVector();
      G4int n = vmsc.size();
      for(G4int i=0; i<n; i++) {
        if((vmsc[i])->GetProcessName() == currentName && (vmsc[i])->Particle() == p) {
          currentLambda = (vmsc[i])->LambdaTable();
	  break;
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::FindEmModel(const G4ParticleDefinition* p,
                                 const G4String& processName,
				       G4double kinEnergy, size_t idxRegion)
{
  // Search for the process
  currentName = processName;

  G4LossTableManager* manager = G4LossTableManager::Instance();
  const std::vector<G4VEnergyLossProcess*> vel = manager->GetEnergyLossProcessVector();
  G4int n = vel.size();
  for(G4int i=0; i<n; i++) {
    if((vel[i])->GetProcessName() == currentName && (vel[i])->Particle() == p) {
      const G4ParticleDefinition* bp = (vel[i])->BaseParticle();
      if(bp) {
        massRatio = bp->GetPDGMass()/p->GetPDGMass();
	baseParticle = bp;
      }
      currentModel = (vel[i])->SelectModelForMaterial(kinEnergy*massRatio, idxRegion);
      break;
    }
  }
  if(!currentModel) {
    const std::vector<G4VEmProcess*> vem = manager->GetEmProcessVector();
    G4int n = vem.size();
    for(G4int i=0; i<n; i++) {
      if((vem[i])->GetProcessName() == currentName && (vem[i])->Particle() == p) {
        currentModel = (vem[i])->SelectModelForMaterial(kinEnergy, idxRegion);
        break;
      }
    }
  }
  if(!currentLambda) {
    const std::vector<G4VMultipleScattering*> vmsc = manager->GetMultipleScatteringVector();
    G4int n = vmsc.size();
    for(G4int i=0; i<n; i++) {
      if((vmsc[i])->GetProcessName() == currentName && (vmsc[i])->Particle() == p) {
        currentModel = (vmsc[i])->SelectModelForMaterial(kinEnergy, idxRegion);
        break;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4VEnergyLossProcess* G4EmCalculator::FindEnergyLossProcess(
                      const G4ParticleDefinition* p)
{
  const G4VEnergyLossProcess* elp = 0;
  G4LossTableManager* manager = G4LossTableManager::Instance();
  const std::vector<G4VEnergyLossProcess*> vel = manager->GetEnergyLossProcessVector();
  G4int n = vel.size();
  for(G4int i=0; i<n; i++) {
    if((vel[i])->Particle() == p) {
      elp = vel[i];
      break;
    }
  }
  return elp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

