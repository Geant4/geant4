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
// GEANT4 Class file
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
// 23-12-02 V.Ivanchenko change interface in order to move
//                           to cut per region
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
// 24-01-03 Make models region aware (V.Ivanchenko)
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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4EmModelManager.hh"
#include "G4LossTableManager.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"
#include "G4DataVector.hh"
#include "G4PhysicsVector.hh"
#include "G4VParticleChange.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmModelManager::G4EmModelManager():
  nEmModels(0),
  nmax(4),
  orderIsChanged(false),
  minSubRange(0.1),
  particle(0)
{
  verboseLevel = 0;
  for(G4int i = 0; i<nmax; i++) {
    emModels[i] = 0;
    order[i] = 0;
    upperEkin[i] = 0.0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmModelManager::~G4EmModelManager()
{
  Clear();
  for(G4int i = 0; i<nmax; i++) {
    if(emModels[i]) delete emModels[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::Clear()
{
  if(0 < verboseLevel) {
    G4cout << "G4EmModelManager::Clear()" << G4endl;
  }

  theCuts.clear();
  theSubCuts.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4DataVector* G4EmModelManager::Initialise(const G4ParticleDefinition* p,
                                                 const G4ParticleDefinition* sp,
						       G4double theMinSubRange,
                                                       G4int val)
{
  // Are models defined?
  if(!nEmModels) {
    G4Exception("G4EmModelManager::Initialise without any model defined");
  }
  particle = p;
  secondaryParticle = sp;
  minSubRange = theMinSubRange;
  verboseLevel = val;

  if(0 < verboseLevel) {
    G4cout << "G4EmModelManager::Initialise() for "
           << p->GetParticleName()
           << G4endl;
  }

  // Ordering
  if(orderIsChanged) {
    G4int oldOrder[5];
    G4VEmModel* oldEmModels[5];
    for(G4int k = 0; k<nmax; k++) {
      oldEmModels[k] = emModels[k];
      oldOrder[k] = order[k];
    }


    for(G4int ik=0; ik<nEmModels; ik++) {

      G4int low = INT_MAX;
      G4int index = 0;
      for(G4int kk=0; kk<nEmModels; kk++) {
        if(oldEmModels[kk]) {
          if(oldOrder[kk] < low) {
            low = oldOrder[kk];
            index = kk;
	  }
	}
      }
      emModels[ik] = oldEmModels[index];
      order[ik] = ik;
      oldEmModels[index] = 0;
    }
    orderIsChanged = false;
  }

  G4DataVector eLow;
  eLow.clear();

  G4int n = 0;

  for(G4int j=0; j<nEmModels; j++) {

    G4double ep = 0.0;
    if(0 < n) ep = upperEkin[n-1];
    G4VEmModel* model = emModels[j];
    G4bool accepted = false;

    if(model->IsInCharge(particle)) {

      G4double tmin = model->LowEnergyLimit(particle);
      G4double tmax = model->HighEnergyLimit(particle);

      if(1 < verboseLevel) {
          G4cout << "New model for tmin(MeV)= " << tmin/MeV
                 << "; tmax(MeV)= " << tmax/MeV
                 << "; tlast(MeV)= " << ep/MeV
                 << G4endl;
      }

      if(tmax > tmin) {

        if(n == 0 || tmax > upperEkin[n-1]) {

          // First model or next model for more high energy range;
          upperEkin[n] = (tmax);
          if(0 < n) tmin = G4std::max(tmin, upperEkin[n-1]);
          eLow.push_back(tmin);
          n++;
          accepted = true;

        } else {

          G4cout << "The model number #" << j
                 << " has no active range "
                 << "; tmax(MeV)= " << tmax/MeV
                 << "; tlast(MeV)= " << ep/MeV
                 << G4endl;
	}

      } else {

          G4cout << "The model number #" << j
                 << " has no active range "
                 << "; tmin(MeV)= " << tmin/MeV
                 << "; tmax(MeV)= " << tmax/MeV
                 << G4endl;
      }
    }
    if(!accepted) {
      upperEkin[n] = ep;
      eLow.push_back(ep);
      n++;
    }
  }

  // Access to materials and build cuts
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  for(size_t i=0; i<numOfCouples; i++) {

    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
    const G4Material* material = couple->GetMaterial();

    if(1 < verboseLevel) {
      G4cout << "G4EmModelManager::Initialise() for "
             << material->GetName() << G4endl;
    }

    G4double cut = 0.0;
    G4double subcut = 0.0;
    if(secondaryParticle) {
      size_t idx = 1;
      if( secondaryParticle == G4Gamma::Gamma() ) idx = 0;
      cut = (*theCoupleTable->GetEnergyCutsVector(idx))[i];
      subcut = minSubRange*cut;
    }

    for(G4int j=0; j<nEmModels; j++) {

      G4VEmModel* model = emModels[j];
      if(upperEkin[j] > eLow[j]) {

        G4double tcutmin = model->MinEnergyCut(particle, couple);

        if(1 < verboseLevel) {
            G4cout << "The model # " << j
                   << "; tcutmin(MeV)= " << tcutmin/MeV
                   << G4endl;
        }

        cut = G4std::max(cut, tcutmin);
	G4double x = G4std::max(cut*minSubRange, tcutmin);
        subcut = G4std::max(subcut, x);
      }
    }
    theCuts.push_back(cut);
    theSubCuts.push_back(subcut);
  }

  for(G4int jj=0; jj<nEmModels; jj++) {
    emModels[jj]->Initialise(particle, theCuts);
  }

  if(1 < verboseLevel) {
    G4cout << "G4EmModelManager is initialised "
           << G4endl;
  }

  return &theCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::AddEmModel(G4VEmModel* p, G4int num, const G4Region*)
{
  if(nEmModels) {
    for(G4int i=0; i<nEmModels; i++) {
      if(num < order[i]) orderIsChanged = true;
    }
  }

  if(nEmModels == nmax) {

    G4cout << "G4EmModelManager::AddEmModel WARNING: cannot accept model #"
           << nEmModels << " - the list is closed"
           << G4endl;

  } else {

    emModels[nEmModels] = p;
    order[nEmModels] = num;
    currentModel = p;
    nEmModels++;

  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::FillDEDXVector(G4PhysicsVector* aVector,
                                const G4MaterialCutsCouple* couple)
{

  // vectors to provide continues dE/dx
  G4DataVector factor;
  G4DataVector dedxLow;
  G4DataVector dedxHigh;

  G4double e;

  const G4Material* material = couple->GetMaterial();
  size_t i = couple->GetIndex();
  G4double cut = theCuts[i];

  if(0 < verboseLevel) {
    G4cout << "G4EmModelManager::FillDEDXVector() for "
           << material->GetName()
	   << "  Ecut(MeV)= " << cut/MeV
           << G4endl;
  }

  factor.resize(nEmModels);
  dedxLow.resize(nEmModels);
  dedxHigh.resize(nEmModels);

  if(0 < verboseLevel) {
      G4cout << "There are " << nEmModels << " models for "
             << material->GetName() << G4endl;
  }

  // calculate factors to provide continuity of energy loss
  factor[0] = 1.0;
  G4int j;
  G4int totBinsLoss = aVector->GetVectorLength();

  dedxLow[0]  = 0.0;

  e = upperEkin[0];
  dedxHigh[0] = emModels[0]->ComputeDEDX(material,particle,e,cut);

  if(nEmModels > 1) {
    for(j=1; j<nEmModels; j++) {

      e = upperEkin[j-1];
      dedxLow[j] = emModels[j]->ComputeDEDX(material,particle,e,cut);
      e = upperEkin[j];
      dedxHigh[j] = emModels[j]->ComputeDEDX(material,particle,e,cut);
    }

    for(j=1; j<nEmModels; j++) {
      if(dedxLow[j] > 0.0) factor[j] = (dedxHigh[j-1]/dedxLow[j] - 1.0);
      else                 factor[j] = 0.0;
    }

    if(0 < verboseLevel) {
      G4cout << "Loop over " << totBinsLoss << " bins start " << G4endl;
    }
  }

  // Calculate energy losses vector
  for(j=0; j<totBinsLoss; j++) {

    G4double e = aVector->GetLowEdgeEnergy(j);
    G4double fac = 1.0;

      // Choose a model of energy losses
    G4int k = 0;
    if(nEmModels > 1) {
      if (e >= upperEkin[0]) {
        for(k=1; k<nEmModels; k++) {
          fac *= (1.0 + factor[k]*upperEkin[k-1]/e);
          if(e <= upperEkin[k] || k == nEmModels-1) break;
        }
      }
    }

    G4double dedx = emModels[k]->ComputeDEDX(material,particle,e,cut)*fac;

    if(dedx < 0.0) dedx = 0.0;
    if(0 < verboseLevel) {
        G4cout << "Material= " << material->GetName()
               << "   E(MeV)= " << e/MeV
               << "  dEdx(MeV/mm)= " << dedx*mm/MeV
               << G4endl;
    }
    aVector->PutValue(j, dedx);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::FillLambdaVector(G4PhysicsVector* aVector,
                                  const G4MaterialCutsCouple* couple)
{

  if(0 < verboseLevel) {
    G4cout << "G4EmModelManager::FillLambdaVector() for particle "
           << particle->GetParticleName() << G4endl;
  }


  // vectors to provide continues dE/dx
  G4DataVector factor;
  G4DataVector sigmaLow;
  G4DataVector sigmaHigh;

  G4double e;

  const G4Material* material = couple->GetMaterial();
  size_t i = couple->GetIndex();
  G4double cut = theCuts[i];
  factor.resize(nEmModels);
  sigmaLow.resize(nEmModels);
  sigmaHigh.resize(nEmModels);

  if(0 < verboseLevel) {
      G4cout << "There are " << nEmModels << " models for "
             << material->GetName() << G4endl;
  }

  // calculate factors to provide continuity of energy loss
  factor[0] = 1.0;
  G4int j;
  G4int totBinsLambda = aVector->GetVectorLength();

  sigmaLow[0]  = 0.0;

  e = upperEkin[0];

  if(1 < verboseLevel) {
      G4cout << "### For material " << material->GetName()
             << "  " << nEmModels
             << " models"
             << " Ecut(MeV)= " << cut/MeV
             << " Emax(MeV)= " << e/MeV
             << " nbins= "  << totBinsLambda
             << G4endl;
  }

  sigmaHigh[0] = emModels[0]->CrossSection(material,particle,e,cut,e);

  if(nEmModels > 1) {

    for(j=1; j<nEmModels; j++) {

      e  = upperEkin[j-1];
      sigmaLow[j] = emModels[j]->CrossSection(material,particle,e,cut,e);
      e  = upperEkin[j];
      sigmaHigh[j] = emModels[j]->CrossSection(material,particle,e,cut,e);
    }
    for(j=1; j<nEmModels; j++) {
      if(sigmaLow[j] > 0.0) factor[j] = (sigmaHigh[j-1]/sigmaLow[j] - 1.0);
      else                  factor[j] = 0.0;
    }
  }

  // Calculate lambda vector
  for(j=0; j<totBinsLambda; j++) {

    e = aVector->GetLowEdgeEnergy(j);

    // Choose a model of energy losses
    G4int k = 0;
    G4double fac = 1.0;
    if(nEmModels > 1) {
      if(e >= upperEkin[0]) {
        for(k=1; k<nEmModels; k++) {
          fac *= (1.0 + factor[k]*upperEkin[k-1]/e);
          if(e <= upperEkin[k] || k == nEmModels-1) break;
        }
      }
    }

    // Cross section interpolation should start from zero
    G4double cross = 0.0;
    if(j > 0) {
      cross = emModels[k]->CrossSection(material,particle,e,cut,e)*fac;
    }

    if(1 < verboseLevel) {
      G4cout << "BuildLambdaTable: " << j << ".   e(MeV)= " << e/MeV
               << "  cross(1/mm)= " << cross*mm
               << " fac= " << fac
               << G4endl;
    }
    if(cross <= 0.0) cross = 0.0;
    //    if(cross <= 0.0) cross = DBL_MAX;
    //    else             cross = 1.0/cross;

    aVector->PutValue(j, cross);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::FillSubLambdaVector(G4PhysicsVector* aVector,
                                     const G4MaterialCutsCouple* couple)
{
  if(0 < verboseLevel) {
    G4cout << "G4EmModelManager::BuildLambdaSubTable() for particle "
           << particle->GetParticleName() << G4endl;
  }


  // vectors to provide continues dE/dx
  G4DataVector factor;
  G4DataVector sigmaLow;
  G4DataVector sigmaHigh;

  G4double e;

  const G4Material* material = couple->GetMaterial();
  size_t i = couple->GetIndex();
  G4double cut = theCuts[i];
  G4double subcut = theSubCuts[i];
  factor.resize(nEmModels);
  sigmaLow.resize(nEmModels);
  sigmaHigh.resize(nEmModels);

  if(0 < verboseLevel) {
      G4cout << "There are " << nEmModels << " models for "
             << material->GetName() << G4endl;
  }

  // calculate factors to provide continuity of energy loss
  factor[0] = 1.0;
  G4int j;
  G4int totBinsLambda = aVector->GetVectorLength();

  sigmaLow[0]  = 0.0;

  e = upperEkin[0];


  if(1 < verboseLevel) {
      G4cout << "### For material " << material->GetName()
             << " are available " << nEmModels
             << " models"
             << " Ecut(MeV)= " << cut/MeV
             << " nbins= "  << totBinsLambda
             << G4endl;
  }

  sigmaHigh[0] = emModels[0]->CrossSection(material,particle,e,subcut,cut);

  if(nEmModels > 1) {

    for(j=1; j<nEmModels; j++) {

      e  = upperEkin[j-1];
      sigmaLow[j] = emModels[j]->CrossSection(material,particle,e,subcut,cut);
      e  = upperEkin[j];
      sigmaHigh[j] = emModels[j]->CrossSection(material,particle,e,subcut,cut);
    }
    for(j=1; j<nEmModels; j++) {
      if(sigmaLow[j] > 0.0) factor[j] = (sigmaHigh[j-1]/sigmaLow[j] - 1.0);
      else                  factor[j] = 0.0;
    }
  }

  // Calculate energy losses vector
  for(j=0; j<totBinsLambda; j++) {

    e = aVector->GetLowEdgeEnergy(j);

    // Choose a model of energy losses
    G4int k = 0;
    G4double fac = 1.0;
    if(nEmModels > 1) {
      if(e >= upperEkin[0]) {
        for(k=1; k<nEmModels; k++) {
          fac *= (1.0 + factor[k]*upperEkin[k-1]/e);
          if(e <= upperEkin[k] || k == nEmModels-1) break;
        }
      }
    }

    G4double cross = 0.0;

    // Cross section interpolation should start from zero
    if (j > 0) {
      cross = emModels[k]->CrossSection(material,particle,e,subcut,cut)*fac;
    }

    if(1 < verboseLevel) {
        G4cout << "BuildLambdaTable: e(MeV)= " << e/MeV
               << "  cross(1/mm)= " << cross*mm
               << " fac= " << fac
               << G4endl;
    }
    if(cross <= 0.0) cross = 0.0;

    aVector->PutValue(j, cross);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
