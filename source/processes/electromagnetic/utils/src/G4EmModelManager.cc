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
// 13-02-03 The set of models is defined for region (V.Ivanchenko)
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
#include "G4DataVector.hh"
#include "G4PhysicsVector.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4RegionModels::G4RegionModels(G4int nMod, G4std::vector<G4int>& list, G4DataVector& lowE)
{
  nModelsForRegion      = nMod;
  theListOfModelIndexes = new G4int [nModelsForRegion];
  lowKineticEnergy      = new G4double [nModelsForRegion];
  for (G4int i=0; i<nModelsForRegion; i++) {
    theListOfModelIndexes[i] = list[i];
    lowKineticEnergy[i] = lowE[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4RegionModels::~G4RegionModels()
{
  delete [] theListOfModelIndexes;
  delete [] lowKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmModelManager::G4EmModelManager():
  nEmModels(0),
  nRegions(0),
  nCouples(0),
  minSubRange(0.1),
  particle(0),
  verboseLevel(0)
{
  models.clear();
  flucModels.clear();
  regions.clear();
  orderOfModels.clear();
  upperEkin.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmModelManager::~G4EmModelManager()
{
  G4int i,j;
  Clear();
  for(i = 0; i<nEmModels; i++) {
    orderOfModels[i] = 1;
  }
  for(i = 0; i<nEmModels; i++) {
    if (orderOfModels[i]) {
      orderOfModels[i] = 0;
      for(j = i+1; j<nEmModels; j++) {
        if(models[i] == models[j]) orderOfModels[j] = 0;
      }
      delete models[i];
    }
  }
  for(i = 0; i<nEmModels; i++) {
    orderOfModels[i] = 1;
  }
  for(i = 0; i<nEmModels; i++) {
    if (orderOfModels[i]) {
      orderOfModels[i] = 0;
      for(j = i+1; j<nEmModels; j++) {
        if(flucModels[i] == flucModels[j]) orderOfModels[j] = 0;
      }
      delete flucModels[i];
    }
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
  upperEkin.clear();
  idxOfRegionModels.clear();
  setOfRegionModels.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::AddEmModel(G4int num, G4VEmModel* p,
                                  G4VEmFluctuationModel* fm, const G4Region* r)
{
  if(!p) {
    G4cout << "G4EmModelManager::AddEmModel WARNING: no model defined." << G4endl;
    return;
  }
  models.push_back(p);
  flucModels.push_back(fm);
  regions.push_back(r);
  orderOfModels.push_back(num);
  if (nEmModels) {
    G4int idx = nEmModels;
    do {idx--;} while (idx && num < orderOfModels[idx]);
    if (num >= orderOfModels[idx] && num <= orderOfModels[idx+1]) idx++;
    if (idx < nEmModels) {
      models[nEmModels] = models[idx];
      flucModels[nEmModels] = flucModels[idx];
      regions[nEmModels] = regions[idx];
      orderOfModels[nEmModels] = orderOfModels[idx];
      models[idx] = p;
      flucModels[idx] = fm;
      regions[idx] = r;
      orderOfModels[idx] = num;
    }
  }
  nEmModels++;
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
    G4cout << "### G4EmModelManager::Initialise() for "
           << p->GetParticleName()
           << G4endl;
  }
  Clear();

  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  const G4Region* world = regionStore->GetRegion("DefaultRegionForTheWorld", false);

  // Identify the list of regions with different set of models
  nRegions = 1;
  G4std::vector<const G4Region*> set;
  set.push_back(world);

  for (G4int ii=0; ii<nEmModels; ii++) {
    const G4Region* r = regions[ii];
    if ( r && r != world) {
      G4bool newRegion = true;
      if (nRegions>1) {
        for (G4int j=1; j<nRegions; j++) {
	  if ( r == set[j] ) newRegion = false;
        }
      }
      if (newRegion) {
        set.push_back(r);
	nRegions++;
      }
    }
  }

  setOfRegionModels.clear();
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  idxOfRegionModels.resize(numOfCouples);
  upperEkin.resize(numOfCouples);

  // Order models for regions
  for (G4int reg=0; reg<nRegions; reg++) {

    const G4Region* region = set[reg];

    G4int n = 0;

    G4std::vector<G4int>    modelAtRegion;
    G4DataVector            eLow;
    G4DataVector            eHigh;
    modelAtRegion.clear();
    eLow.clear();
    eHigh.clear();

    for (G4int ii=0; ii<nEmModels; ii++) {

      G4VEmModel* model = models[ii];
      if ( (model->IsInCharge(particle)) &&
           (0 == regions[ii] || region == regions[ii]) )
      {

        G4double tmin = model->LowEnergyLimit(particle);
        G4double tmax = model->HighEnergyLimit(particle);

        if(0 < verboseLevel) {
          G4cout << "Model # " << ii << " for region <"
	         << region->GetName() << ">  "
	         << " tmin(MeV)= " << tmin/MeV
                 << "; tmax(MeV)= " << tmax/MeV
                 << G4endl;
        }
        if (n) tmin = G4std::max(tmin, eHigh[n]);

	if (tmin < tmax) {
          modelAtRegion.push_back(ii);
 	  eLow.push_back(tmin);
	  eHigh.push_back(tmax);
	  upperEkin[ii] = tmax;
	  n++;
	}
	eLow[0] = 0.0;
      }
    }

    if(0 < verboseLevel) {
      G4cout << "New G4RegionModels set with " << n << " models for region <"
	     << region->GetName() << ">  " << G4endl;
    }
    G4RegionModels* rm = new G4RegionModels(n, modelAtRegion, eLow);
    setOfRegionModels.push_back(rm);
  }

  // Access to materials and build cuts

  for(size_t i=0; i<numOfCouples; i++) {

    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
    const G4Material* material = couple->GetMaterial();
    const G4ProductionCuts* pcuts = couple->GetProductionCuts();
    G4int reg = nRegions;
    do {reg--;} while (reg && pcuts != (set[reg]->GetProductionCuts()));
    idxOfRegionModels[i] = reg;

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

    G4int nm = setOfRegionModels[reg]->NumberOfModels();
    for(G4int j=0; j<nm; j++) {

      G4VEmModel* model = models[setOfRegionModels[reg]->ModelIndex(j)];

      G4double tcutmin = model->MinEnergyCut(particle, couple);

      cut = G4std::max(cut, tcutmin);
      G4double x = G4std::max(cut*minSubRange, tcutmin);
      subcut = G4std::max(subcut, x);
      if(0 < verboseLevel) {
            G4cout << "The model # " << j
                   << "; tcutmin(MeV)= " << tcutmin/MeV
                   << "; tcut(MeV)= " << cut/MeV
                    << G4endl;
      }
    }
    theCuts.push_back(cut);
    theSubCuts.push_back(subcut);
  }

  for(G4int jj=0; jj<nEmModels; jj++) {
    models[jj]->Initialise(particle, theCuts);
    if(flucModels[jj]) flucModels[jj]->Initialise(particle);
  }


  if(0 < verboseLevel) {
    G4cout << "G4EmModelManager is initialised "
           << G4endl;
  }

  return &theCuts;
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

  G4int reg  = idxOfRegionModels[i];
  const G4RegionModels* regModels = setOfRegionModels[reg];
  G4int nmod = regModels->NumberOfModels();
  factor.resize(nmod);
  dedxLow.resize(nmod);
  dedxHigh.resize(nmod);


  if(0 < verboseLevel) {
      G4cout << "There are " << nmod << " models for "
             << material->GetName()
	     << " at the region #" << reg
	     << G4endl;
  }


  // calculate factors to provide continuity of energy loss
  factor[0] = 1.0;
  G4int j;

  G4int totBinsLoss = aVector->GetVectorLength();

  dedxLow[0]  = 0.0;

  e = upperEkin[regModels->ModelIndex(0)];
  dedxHigh[0] = models[regModels->ModelIndex(0)]->ComputeDEDX(material,particle,e,cut);

  if(nmod > 1) {
    for(j=1; j<nmod; j++) {

      e = upperEkin[regModels->ModelIndex(j-1)];
      dedxLow[j] = models[regModels->ModelIndex(j)]->ComputeDEDX(material,particle,e,cut);
      e = upperEkin[regModels->ModelIndex(j)];
      dedxHigh[j] = models[regModels->ModelIndex(j)]->ComputeDEDX(material,particle,e,cut);
    }

    for(j=1; j<nmod; j++) {
      if(dedxLow[j] > 0.0) factor[j] = (dedxHigh[j-1]/dedxLow[j] - 1.0);
      else                 factor[j] = 0.0;
    }

    if(1 < verboseLevel) {
      G4cout << "Loop over " << totBinsLoss << " bins start " << G4endl;
    }
  }

  // Calculate energy losses vector
  for(j=0; j<totBinsLoss; j++) {

    G4double e = aVector->GetLowEdgeEnergy(j);
    G4double fac = 1.0;

      // Choose a model of energy losses
    G4int k = 0;
    if (nmod > 1 && e > upperEkin[regModels->ModelIndex(0)]) {
      do {
        k++;
        fac *= (1.0 + factor[k]*upperEkin[regModels->ModelIndex(k-1)]/e);
      } while (k<nmod-1 && e < upperEkin[regModels->ModelIndex(k)] );
    }

    G4double dedx = models[regModels->ModelIndex(k)]->ComputeDEDX(material,particle,e,cut)*fac;

    if(dedx < 0.0) dedx = 0.0;
    if(1 < verboseLevel) {
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

  G4int reg  = idxOfRegionModels[i];
  const G4RegionModels* regModels = setOfRegionModels[reg];
  G4int nmod = regModels->NumberOfModels();
  factor.resize(nmod);
  sigmaLow.resize(nmod);
  sigmaHigh.resize(nmod);

  if(0 < verboseLevel) {
      G4cout << "There are " << nmod << " models for "
             << material->GetName() << G4endl;
  }

  // calculate factors to provide continuity of energy loss
  factor[0] = 1.0;
  G4int j;
  G4int totBinsLambda = aVector->GetVectorLength();

  sigmaLow[0]  = 0.0;

  e = upperEkin[regModels->ModelIndex(0)];

  if(1 < verboseLevel) {
      G4cout << "### For material " << material->GetName()
             << "  " << nmod
             << " models"
             << " Ecut(MeV)= " << cut/MeV
             << " Emax(MeV)= " << e/MeV
             << " nbins= "  << totBinsLambda
             << G4endl;
  }

  sigmaHigh[0] = models[regModels->ModelIndex(0)]->CrossSection(material,particle,e,cut,e);

  if(nmod > 1) {

    for(j=1; j<nmod; j++) {

      e  = upperEkin[regModels->ModelIndex(j-1)];
      sigmaLow[j] = models[regModels->ModelIndex(j)]->CrossSection(material,particle,e,cut,e);
      e  = upperEkin[regModels->ModelIndex(j)];
      sigmaHigh[j] = models[regModels->ModelIndex(j)]->CrossSection(material,particle,e,cut,e);
    }
    for(j=1; j<nmod; j++) {
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
    if (nmod > 1 && e > upperEkin[regModels->ModelIndex(0)]) {
      do {
        k++;
        fac *= (1.0 + factor[k]*upperEkin[regModels->ModelIndex(k-1)]/e);
      } while (k<nmod-1 && e < upperEkin[regModels->ModelIndex(k)] );
    }

    // Cross section interpolation should start from zero
    G4double cross = 0.0;
    if(j > 0) {
      cross = models[regModels->ModelIndex(k)]->CrossSection(material,particle,e,cut,e)*fac;
    }

    if(1 < verboseLevel) {
      G4cout << "BuildLambdaTable: " << j << ".   e(MeV)= " << e/MeV
               << "  cross(1/mm)= " << cross*mm
               << " fac= " << fac
               << G4endl;
    }
    if(cross <= 0.0) cross = 0.0;

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

  G4int reg  = idxOfRegionModels[i];
  const G4RegionModels* regModels = setOfRegionModels[reg];
  G4int nmod = regModels->NumberOfModels();
  factor.resize(nmod);
  sigmaLow.resize(nmod);
  sigmaHigh.resize(nmod);

  if(0 < verboseLevel) {
      G4cout << "There are " << nmod << " models for "
             << material->GetName() << G4endl;
  }

  // calculate factors to provide continuity of energy loss
  factor[0] = 1.0;
  G4int j;
  G4int totBinsLambda = aVector->GetVectorLength();

  sigmaLow[0]  = 0.0;

  e = upperEkin[regModels->ModelIndex(0)];


  if(1 < verboseLevel) {
      G4cout << "### For material " << material->GetName()
             << " are available " << nmod
             << " models"
             << " Ecut(MeV)= " << cut/MeV
             << " nbins= "  << totBinsLambda
             << G4endl;
  }

  sigmaHigh[0] = models[regModels->ModelIndex(0)]->CrossSection(material,particle,e,subcut,cut);

  if(nmod > 1) {

    for(j=1; j<nmod; j++) {

      e  = upperEkin[regModels->ModelIndex(j-1)];
      sigmaLow[j] = models[regModels->ModelIndex(j)]->CrossSection(material,particle,e,subcut,cut);
      e  = upperEkin[regModels->ModelIndex(j)];
      sigmaHigh[j] = models[regModels->ModelIndex(j)]->CrossSection(material,particle,e,subcut,cut);
    }
    for(j=1; j<nmod; j++) {
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
    if (nmod > 1 && e > upperEkin[regModels->ModelIndex(0)]) {
      do {
        k++;
        fac *= (1.0 + factor[k]*upperEkin[regModels->ModelIndex(k-1)]/e);
      } while (k<nmod-1 && e < upperEkin[regModels->ModelIndex(k)] );
    }

    G4double cross = 0.0;

    // Cross section interpolation should start from zero
    if (j > 0) {
      cross = models[regModels->ModelIndex(k)]->CrossSection(material,particle,e,subcut,cut)*fac;
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
