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
// $Id: G4EmCalculator.cc,v 1.1 2004-06-30 14:36:51 vnivanch Exp $
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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmCalculator::G4EmCalculator()
{
  theManager = G4LossTableManager::Instance();
  nLocalMaterials    = 0;
  currentCoupleIndex = 0;
  currentCouple      = 0;
  currentMaterial    = 0;
  currentParticle    = 0;
  currentLambda      = 0;
  chargeSquare       = 1.0;
  massRatio          = 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmCalculator::~G4EmCalculator()
{
  for (G4int i=0; i<nLocalMaterials; i++) {
    delete localCouples[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetDEDX(const G4ParticleDefinition* p, const G4Material* mat,
                                       G4double kinEnergy)
{
  G4double res = 0.0;
  FindMaterial(mat, p);
  if(currentCouple) res = theManager->GetDEDX(p, kinEnergy, currentCouple);
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetRange(const G4ParticleDefinition* p, const G4Material* mat,
                                        G4double kinEnergy)
{
  G4double res = 0.0;
  FindMaterial(mat, p);
  if(currentCouple) res = theManager->GetRange(p, kinEnergy, currentCouple);
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetKinEnergy(const G4ParticleDefinition* p, const G4Material* mat,
                                        G4double range)
{
  G4double res = 0.0;
  FindMaterial(mat, p);
  if(currentCouple) res = theManager->GetEnergy(p, range, currentCouple);
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetCrossSection(const G4ParticleDefinition* p, const G4Material* mat,
                                         const G4String& processName, G4double kinEnergy)
{
  G4double res = 0.0;
  FindMaterial(mat, p);
  if(currentCouple) {
    FindLambdaTable(p, processName);
    if(currentLambda) {
      G4bool b;
      res = (((*currentLambda)[currentCoupleIndex])->GetValue(kinEnergy*massRatio,b))*chargeSquare;
      res /= currentMaterial->GetTotNbOfAtomsPerVolume();
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::GetMinFreePath(const G4ParticleDefinition* p, const G4Material* mat,
                                        const G4String& processName, G4double kinEnergy)
{
  G4double res = DBL_MAX;
  G4double x = GetCrossSection(p, mat, processName, kinEnergy);
  if(x > 0.0) res = 1.0/(x*(currentMaterial->GetTotNbOfAtomsPerVolume()));
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeDEDX(const G4ParticleDefinition* p, const G4Material* mat,
                                     const G4String& processName,
                                           G4double kinEnergy, G4double cut, size_t idxRegion)
{
  G4double res = 0.0;
  FindCouple(mat, cut);
  FindEmModel(p, processName, kinEnergy, idxRegion);
  if(currentModel)
    res = currentModel->ComputeDEDX(currentCouple, baseParticle, kinEnergy*massRatio, currentCut);
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeCrossSection(const G4ParticleDefinition* p, const G4Material* mat,
                                             const G4String& processName,
                                           G4double kinEnergy, G4double cut, size_t idxRegion)
{
  G4double res = 0.0;
  FindCouple(mat, cut);
  FindEmModel(p, processName, kinEnergy, idxRegion);
  if(currentModel) {
    G4double e = kinEnergy*massRatio;
    res = currentModel->CrossSection(currentCouple, baseParticle, e, currentCut, e);
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmCalculator::ComputeMinFreePath(const G4ParticleDefinition* p, const G4Material* mat,
                                            const G4String& processName,
                                           G4double kinEnergy, G4double cut, size_t idxRegion)
{
  G4double mfp = DBL_MAX;
  G4double x = ComputeCrossSection(p, mat, processName, kinEnergy, cut, idxRegion);
  if(x > 0.0) mfp = 1.0/(x*(currentMaterial->GetTotNbOfAtomsPerVolume()));
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmCalculator::FindMaterial(const G4Material* material, const G4ParticleDefinition* p)
{
  if(currentMaterial == material) return;
  baseParticle = p;
  G4double q = p->GetPDGCharge()/eplus;
  chargeSquare = q*q;

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
    currentParticle = p;
    baseParticle    = p;
    currentName     = processName;
    currentLambda   = 0;
    massRatio       = 1.0;
    G4LossTableManager* theManager = G4LossTableManager::Instance();
    const std::vector<G4VEnergyLossProcess*> vel = theManager->GetEnergyLossProcessVector();
    G4int n = vel.size();
    for(G4int i=0; i<n; i++) {
      if((vel[i])->GetProcessName() == currentName) {
        currentLambda = (vel[i])->LambdaTable();
	const G4ParticleDefinition* bp = (vel[i])->BaseParticle();
	if(bp) {
	  massRatio = bp->GetPDGMass()/p->GetPDGMass();
	  baseParticle = bp;
	}
	break;
      }
    }
    if(!currentLambda) {
      const std::vector<G4VEmProcess*> vem = theManager->GetEmProcessVector();
      G4int n = vem.size();
      for(G4int i=0; i<n; i++) {
        if((vem[i])->GetProcessName() == currentName) {
          currentLambda = (vem[i])->LambdaTable();
	  break;
        }
      }
    }
    if(!currentLambda) {
      const std::vector<G4VMultipleScattering*> vmsc = theManager->GetMultipleScatteringVector();
      G4int n = vmsc.size();
      for(G4int i=0; i<n; i++) {
        if((vmsc[i])->GetProcessName() == currentName) {
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
  currentParticle = p;
  baseParticle    = p;
  currentName     = processName;
  currentModel    = 0;
  massRatio       = 1.0;
  G4LossTableManager* theManager = G4LossTableManager::Instance();
  const std::vector<G4VEnergyLossProcess*> vel = theManager->GetEnergyLossProcessVector();
  G4int n = vel.size();
  for(G4int i=0; i<n; i++) {
    if((vel[i])->GetProcessName() == currentName) {
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
    const std::vector<G4VEmProcess*> vem = theManager->GetEmProcessVector();
    G4int n = vem.size();
    for(G4int i=0; i<n; i++) {
      if((vem[i])->GetProcessName() == currentName) {
        currentModel = (vem[i])->SelectModelForMaterial(kinEnergy, idxRegion);
        break;
      }
    }
  }
  if(!currentLambda) {
    const std::vector<G4VMultipleScattering*> vmsc = theManager->GetMultipleScatteringVector();
    G4int n = vmsc.size();
    for(G4int i=0; i<n; i++) {
      if((vmsc[i])->GetProcessName() == currentName) {
        currentModel = (vmsc[i])->SelectModelForMaterial(kinEnergy, idxRegion);
        break;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

