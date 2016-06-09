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
// $Id: G4EmModelManager.cc,v 1.46 2008/10/13 14:56:56 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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
// 06-03-03 Fix in energy intervals for models (V.Ivanchenko)
// 13-04-03 Add startFromNull (V.Ivanchenko)
// 13-05-03 Add calculation of precise range (V.Ivanchenko)
// 16-07-03 Replace G4Material by G4MaterialCutCouple in dE/dx and CrossSection
//          calculation (V.Ivanchenko)
// 21-07-03 Add UpdateEmModel method (V.Ivanchenko)
// 03-11-03 Substitute STL vector for G4RegionModels (V.Ivanchenko)
// 26-01-04 Fix in energy range conditions (V.Ivanchenko)
// 24-03-05 Remove check or IsInCharge (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 18-08-05 Fix cut for e+e- pair production (V.Ivanchenko)
// 29-11-05 Add protection for arithmetic operations with cut=DBL_MAX (V.Ivanchenko)
// 20-01-06 Introduce G4EmTableType and reducing number of methods (VI)
// 13-05-06 Add GetModel by index method (VI)
// 15-03-07 Add maxCutInRange (V.Ivanchenko)
// 12-04-07 Add verbosity at destruction (V.Ivanchenko)
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4RegionModels::G4RegionModels(G4int nMod, std::vector<G4int>& indx, 
			       G4DataVector& lowE, const G4Region* reg)
{
  nModelsForRegion      = nMod;
  theListOfModelIndexes = new G4int [nModelsForRegion];
  lowKineticEnergy      = new G4double [nModelsForRegion+1];
  for (G4int i=0; i<nModelsForRegion; i++) {
    theListOfModelIndexes[i] = indx[i];
    lowKineticEnergy[i] = lowE[i];
  }
  lowKineticEnergy[nModelsForRegion] = lowE[nModelsForRegion];
  theRegion = reg;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4RegionModels::~G4RegionModels()
{
  delete [] theListOfModelIndexes;
  delete [] lowKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsVector.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4RegionStore.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmModelManager::G4EmModelManager():
  nEmModels(0),
  nRegions(0),
  nCouples(0),
  idxOfRegionModels(0),
  setOfRegionModels(0),
  minSubRange(0.1),
  particle(0),
  verboseLevel(0)
{
  models.clear();
  flucModels.clear();
  regions.clear();
  orderOfModels.clear();
  maxCutInRange    = 12.*cm;
  maxSubCutInRange = 0.7*mm;
  theGamma = G4Gamma::Gamma();
  thePositron = G4Positron::Positron();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmModelManager::~G4EmModelManager()
{
  Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::Clear()
{
  if(1 < verboseLevel) {
    G4cout << "G4EmModelManager::Clear()" << G4endl;
  }

  theCuts.clear();
  theSubCuts.clear();
  if(idxOfRegionModels) delete [] idxOfRegionModels;
  if(setOfRegionModels && nRegions) {
    for(G4int i=0; i<nRegions; i++) {
      delete (setOfRegionModels[i]);
    }
    delete [] setOfRegionModels;
  }
  idxOfRegionModels = 0;
  setOfRegionModels = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::AddEmModel(G4int num, G4VEmModel* p,
                                  G4VEmFluctuationModel* fm, const G4Region* r)
{
  if(!p) {
    G4cout << "G4EmModelManager::AddEmModel WARNING: no model defined." 
	   << G4endl;
    return;
  }
  models.push_back(p);
  flucModels.push_back(fm);
  regions.push_back(r);
  orderOfModels.push_back(num);
  p->DefineForRegion(r);
  nEmModels++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::UpdateEmModel(const G4String& nam, 
				     G4double emin, G4double emax)
{
  if (nEmModels) {
    for(G4int i=0; i<nEmModels; i++) {
      if(nam == models[i]->GetName()) {
        models[i]->SetLowEnergyLimit(emin);
        models[i]->SetHighEnergyLimit(emax);
	break;
      }
    }
  }
  G4cout << "G4EmModelManager::UpdateEmModel WARNING: no model <"
         << nam << "> is found out"
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmModel* G4EmModelManager::GetModel(G4int i, G4bool ver)
{
  G4VEmModel* m = 0;
  if(i >= 0 && i < nEmModels) {m = models[i];}
  else if(verboseLevel > 0 && ver) { 
    G4cout << "G4EmModelManager::GetModel WARNING: "
	   << "index " << i << " is wrong Nmodels= "
	   << nEmModels;
    if(particle) G4cout << " for " << particle->GetParticleName(); 
    G4cout<< G4endl;
  }
  return m;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4DataVector* G4EmModelManager::Initialise(const G4ParticleDefinition* p,
                                                 const G4ParticleDefinition* sp,
						 G4double theMinSubRange,
						 G4int val)
{
  verboseLevel = val;
  if(1 < verboseLevel) {
    G4cout << "G4EmModelManager::Initialise() for "
           << p->GetParticleName()
           << G4endl;
  }
  // Are models defined?
  if(!nEmModels) {
    G4Exception("G4EmModelManager::Initialise without any model defined");
  }
  particle = p;
  secondaryParticle = sp;
  minSubRange = theMinSubRange;

  Clear();
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  const G4Region* world = 
    regionStore->GetRegion("DefaultRegionForTheWorld", false);

  // Identify the list of regions with different set of models
  nRegions = 1;
  std::vector<const G4Region*> setr;
  setr.push_back(world);
  G4bool isWorld = false;

  for (G4int ii=0; ii<nEmModels; ii++) {
    const G4Region* r = regions[ii];
    if ( r == 0 || r == world) {
      isWorld = true;
      regions[ii] = world;
    } else {
      G4bool newRegion = true;
      if (nRegions>1) {
        for (G4int j=1; j<nRegions; j++) {
	  if ( r == setr[j] ) newRegion = false;
        }
      }
      if (newRegion) {
        setr.push_back(r);
	nRegions++;
      }
    }
  }

  G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = theCoupleTable->GetTableSize();
  idxOfRegionModels = new G4int[numOfCouples+1];
  idxOfRegionModels[numOfCouples] = 0;
  setOfRegionModels = new G4RegionModels*[nRegions];

  std::vector<G4int>    modelAtRegion(nEmModels);
  std::vector<G4int>    modelOrd(nEmModels);
  G4DataVector          eLow(nEmModels);
  G4DataVector          eHigh(nEmModels);
  G4int nmax = nEmModels;

  // Order models for regions
  for (G4int reg=0; reg<nRegions; reg++) {
    const G4Region* region = setr[reg];
    G4int n = 0;

    if(isWorld || 0 < reg) {

      for (G4int ii=0; ii<nEmModels; ii++) {

	G4VEmModel* model = models[ii];
	if ( region == regions[ii] ) {

	  G4double tmin = model->LowEnergyLimit();
	  G4double tmax = model->HighEnergyLimit();
          G4int  ord    = orderOfModels[ii];
          G4bool push   = true;
          G4bool insert = false;
          G4int  idx    = n;

	  if(1 < verboseLevel) {
	    G4cout << "Model #" << ii 
		   << " <" << model->GetName() << "> for region <";
	    if (region) G4cout << region->GetName();
	    G4cout << ">  "
	         << " tmin(MeV)= " << tmin/MeV
                 << "; tmax(MeV)= " << tmax/MeV
		 << "; order= " << ord
                 << G4endl;
	  }
        
	  if (n == 0) n++;
          else {
	    tmin = std::min(tmin, eHigh[n-1]);
	    if(tmin >= tmax) push = false;
            else {

	      // high energy model
	      if(tmin == eHigh[n-1] && tmax > eHigh[n-1]) n++;
	      else if (tmax > eHigh[n-1]) {

		// compare order of models
		for(G4int k = n-1; k>=0; k--) {
       
		  if (ord >= modelOrd[k]) {
		    tmin = std::max(tmin, eHigh[k]);           
		    if(k < n-1) n = k + 2;
		    break;
		  } else if (tmin > eLow[k]) { 
		    eHigh[k] = tmin;
		    n = k + 2;
		    break;
		  } else if (tmin == eLow[k]) {
		    n = k + 1;
		    break;
		  }
		}
		if(tmin < eLow[0]) n = 1;
                idx = n - 1;

		// low energy model
	      } else {

                tmax = std::max(tmax, eLow[0]);
                insert = true;
                push = false;
		idx = 0;
                if(tmax <= eLow[0]) tmax = eLow[0];
		else {

		  for(G4int k=0; k<n; k++) {
       
		    if (ord >= modelOrd[k]) {
                      if(k == 0) {
			if(tmin < eLow[0]) tmax = eLow[0];
			else insert = false;
                        break; 
		      } else {
			insert = false;
                        break;
		      }
		    } else if(tmax < eHigh[k]) {
                      idx = k;
		      if(k > 0) tmin = eLow[k];
                      eLow[k] = tmax;               
		      break; 
		    } else if(tmax == eHigh[k]) {           
		      insert = false;
                      push = true;
                      idx = k;
		      if(k > 0) tmin = eLow[k];
                      else tmin = std::min(tmin,eLow[0]);
		      break;
		    } else {
		      modelAtRegion[k] = ii;
		      modelOrd[k] = ord;
		      if(k == 0) eLow[idx] = std::min(tmin,eLow[0]);
		    }
		  }
		}
                if(insert && idx < n) n++;
                else insert = false;
	      }
	    }
	  }
          if(n > nmax) {
            nmax = n;
            modelAtRegion.resize(nmax);
            modelOrd.resize(nmax);
            eLow.resize(nmax);
            eHigh.resize(nmax);
	  }
	  if(insert) {
            for(G4int k=n-2; k>=idx; k--) {            
	      modelAtRegion[k+1] = modelAtRegion[k];
	      modelOrd[k+1] = modelOrd[k];
	      eLow[k+1]  = eLow[k];
	      eHigh[k+1] = eHigh[k];
	    }
	  }
	  if (push || insert) {
	    modelAtRegion[idx] = ii;
	    modelOrd[idx] = ord;
	    eLow[idx]  = tmin;
	    eHigh[idx] = tmax;
	  }
	}
      }
    } else {
      n = 1;
      models.push_back(0);
      modelAtRegion.push_back(nEmModels);
      eLow.push_back(0.0);
      eHigh.push_back(DBL_MAX);
    }
    eLow[0] = 0.0;
    if(n >= nmax) eLow.resize(nmax+1);
    eLow[n] = eHigh[n-1];

    if(1 < verboseLevel) {
      G4cout << "New G4RegionModels set with " << n << " models for region <";
      if (region) G4cout << region->GetName();
      G4cout << ">  Elow(MeV)= ";
      for(G4int ii=0; ii<=n; ii++) {G4cout << eLow[ii]/MeV << " ";}
      G4cout << G4endl;
    }
    G4RegionModels* rm = new G4RegionModels(n, modelAtRegion, eLow, region);
    setOfRegionModels[reg] = rm;
  }

  // Access to materials and build cuts

  for(G4int i=0; i<numOfCouples; i++) {

    const G4MaterialCutsCouple* couple = 
      theCoupleTable->GetMaterialCutsCouple(i);
    const G4Material* material = couple->GetMaterial();
    const G4ProductionCuts* pcuts = couple->GetProductionCuts();
 
    G4int reg = nRegions;
    do {reg--;} while (reg>0 && pcuts != (setr[reg]->GetProductionCuts()));
    idxOfRegionModels[i] = reg;

    if(1 < verboseLevel) {
      G4cout << "G4EmModelManager::Initialise() for "
             << material->GetName()
             << " indexOfCouple= " << i
             << " indexOfRegion= " << reg
             << G4endl;
    }

    G4double cut = DBL_MAX;
    G4double subcut = DBL_MAX;
    if(secondaryParticle) {
      size_t idx = 1;
      if( secondaryParticle == theGamma ) idx = 0;
      cut = (*theCoupleTable->GetEnergyCutsVector(idx))[i];
      if( secondaryParticle == thePositron && cut < DBL_MAX )
        cut += (*theCoupleTable->GetEnergyCutsVector(2))[i] + 
	       2.0*electron_mass_c2;

      // compute subcut 
      if( cut < DBL_MAX ) subcut = minSubRange*cut;
      if(pcuts->GetProductionCut(idx) < maxCutInRange) {
	G4double tcutmax = 
	  theCoupleTable->ConvertRangeToEnergy(secondaryParticle,
					       material,maxSubCutInRange);
	if(tcutmax < subcut) subcut = tcutmax;
      }
    }

    G4int nm = setOfRegionModels[reg]->NumberOfModels();
    for(G4int j=0; j<nm; j++) {

      G4VEmModel* model = models[setOfRegionModels[reg]->ModelIndex(j)];

      G4double tcutmin = model->MinEnergyCut(particle, couple);

      if(cut < tcutmin) cut = tcutmin;
      if(subcut < tcutmin) subcut = tcutmin;
      if(1 < verboseLevel) {
            G4cout << "The model # " << j
                   << "; tcutmin(MeV)= " << tcutmin/MeV
                   << "; tcut(MeV)= " << cut/MeV
                   << "; tsubcut(MeV)= " << subcut/MeV
                   << " for " << particle->GetParticleName()
		   << G4endl;
      }
    }
    theCuts.push_back(cut);
    theSubCuts.push_back(subcut);
  }

  for(G4int jj=0; jj<nEmModels; jj++) {
    models[jj]->Initialise(particle, theCuts);
    if(flucModels[jj]) flucModels[jj]->InitialiseMe(particle);
  }


  if(1 < verboseLevel) {
    G4cout << "G4EmModelManager for " << particle->GetParticleName() 
           << " is initialised; nRegions=  " << nRegions
           << G4endl;
  }

  return &theCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::FillDEDXVector(G4PhysicsVector* aVector,
				      const G4MaterialCutsCouple* couple,
                                      G4EmTableType tType)
{

  G4double e;

  size_t i = couple->GetIndex();
  G4double cut = theCuts[i];
  G4double subcut = 0.0;

  if(fTotal == tType) cut = DBL_MAX;
  else if(fSubRestricted == tType) subcut = theSubCuts[i];

  if(1 < verboseLevel) {
    G4cout << "G4EmModelManager::FillDEDXVector() for "
           << couple->GetMaterial()->GetName()
	   << "  cut(MeV)= " << cut
	   << "  subcut(MeV)= " << subcut
	   << "  Type " << tType
	   << "  for " << particle->GetParticleName()
           << G4endl;
  }

  G4int reg  = idxOfRegionModels[i];
  const G4RegionModels* regModels = setOfRegionModels[reg];
  G4int nmod = regModels->NumberOfModels();

  // vectors to provide continues dE/dx
  G4DataVector factor(nmod);
  G4DataVector eLow(nmod+1);
  G4DataVector dedxLow(nmod);
  G4DataVector dedxHigh(nmod);

  if(1 < verboseLevel) {
      G4cout << "There are " << nmod << " models for "
             << couple->GetMaterial()->GetName()
	     << " at the region #" << reg
	     << G4endl;
  }

  // calculate factors to provide continuity of energy loss


  factor[0] = 1.0;
  G4int j;

  G4int totBinsLoss = aVector->GetVectorLength();

  dedxLow[0] = 0.0;
  eLow[0]    = 0.0;

  e = regModels->LowEdgeEnergy(1);
  eLow[1]    = e;
  G4VEmModel* model = models[regModels->ModelIndex(0)]; 
  dedxHigh[0] = 0.0;
  if(model && cut > subcut) {
    dedxHigh[0] = model->ComputeDEDX(couple,particle,e,cut);
    if(subcut > 0.0) {
      dedxHigh[0] -= model->ComputeDEDX(couple,particle,e,subcut);
    }
  }
  if(nmod > 1) {
    for(j=1; j<nmod; j++) {

      e = regModels->LowEdgeEnergy(j);
      eLow[j]   = e;
      G4int idx = regModels->ModelIndex(j); 

      dedxLow[j] = models[idx]->ComputeDEDX(couple,particle,e,cut);
      if(subcut > 0.0) {
	dedxLow[j] -= models[idx]->ComputeDEDX(couple,particle,e,subcut);
      }
      if(subcut == cut) dedxLow[j] = 0.0;

      e = regModels->LowEdgeEnergy(j+1);
      eLow[j+1] = e;
      dedxHigh[j] = models[idx]->ComputeDEDX(couple,particle,e,cut);
      if(subcut > 0.0) {
	dedxHigh[j] -= models[idx]->ComputeDEDX(couple,particle,e,subcut);
      }
      if(subcut == cut) dedxHigh[j] = 0.0;
    }
    if(1 < verboseLevel) {
      G4cout << " model #0"  
	     << "  dedx(" << eLow[0] << ")=  " << dedxLow[0]
	     << "  dedx(" << eLow[1] << ")=  " << dedxHigh[0]
	     << G4endl;
    }

    for(j=1; j<nmod; j++) {
      if(dedxLow[j] > 0.0) {
	factor[j] = (dedxHigh[j-1]/dedxLow[j] - 1.0)*eLow[j];
      } else  factor[j] = 0.0;
      if(1 < verboseLevel) {
        G4cout << " model #" << j 
	       << "  dedx(" << eLow[j] << ")=  " << dedxLow[j]
	       << "  dedx(" << eLow[j+1] << ")=  " << dedxHigh[j]
	       << "  factor= " << factor[j]/eLow[j]
	       << G4endl;
      }
    }

    if(2 < verboseLevel) {
      G4cout << "Loop over " << totBinsLoss << " bins start " << G4endl;
    }
  }

  // Calculate energy losses vector
  for(j=0; j<totBinsLoss; j++) {

    G4double e = aVector->GetLowEdgeEnergy(j);
    G4double fac = 1.0;

      // Choose a model of energy losses
    G4int k = 0;
    if (nmod > 1 && e > eLow[1]) {
      do {
        k++;
        fac *= (1.0 + factor[k]/e);
      } while (k+1 < nmod && e > eLow[k+1]);
    }

    model = models[regModels->ModelIndex(k)];
    G4double dedx = 0.0;
    G4double dedx0 = 0.0;
    if(model && cut > subcut) {
      dedx = model->ComputeDEDX(couple,particle,e,cut); 
      dedx0 = dedx;
      if(subcut > 0.0) dedx -= model->ComputeDEDX(couple,particle,e,subcut); 
      dedx *= fac;
    }

    if(dedx < 0.0) dedx = 0.0;
    if(2 < verboseLevel) {
        G4cout << "Material= " << couple->GetMaterial()->GetName()
               << "   E(MeV)= " << e/MeV
               << "  dEdx(MeV/mm)= " << dedx*mm/MeV
               << "  dEdx0(MeV/mm)= " << dedx0*mm/MeV
               << "  fac= " << fac
               << G4endl;
    }
    aVector->PutValue(j, dedx);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::FillLambdaVector(G4PhysicsVector* aVector,
					const G4MaterialCutsCouple* couple,
				        G4bool startFromNull,
					G4EmTableType tType)
{
  G4double e;

  size_t i = couple->GetIndex();
  G4double cut  = theCuts[i];
  G4double tmax = DBL_MAX;
  if (fSubRestricted == tType) {
    tmax = cut;
    cut  = theSubCuts[i];
  }

  if(1 < verboseLevel) {
    G4cout << "G4EmModelManager::FillLambdaVector() for particle "
           << particle->GetParticleName()
           << " in " << couple->GetMaterial()->GetName()
	   << "  Ecut(MeV)= " << cut
	   << "  Emax(MeV)= " << tmax
	   << "  Type " << tType   
           << G4endl;
  }

  G4int reg  = idxOfRegionModels[i];
  const G4RegionModels* regModels = setOfRegionModels[reg];
  G4int nmod = regModels->NumberOfModels();

  // vectors to provide continues dE/dx
  G4DataVector factor(nmod);
  G4DataVector eLow(nmod+1);
  G4DataVector sigmaLow(nmod);
  G4DataVector sigmaHigh(nmod);

  if(2 < verboseLevel) {
      G4cout << "There are " << nmod << " models for "
             << couple->GetMaterial()->GetName() << G4endl;
  }

  // calculate factors to provide continuity of energy loss
  factor[0] = 1.0;
  G4int j;
  G4int totBinsLambda = aVector->GetVectorLength();

  sigmaLow[0] = 0.0;
  eLow[0]     = 0.0;

  e = regModels->LowEdgeEnergy(1);
  eLow[1]     = e;
  G4VEmModel* model = models[regModels->ModelIndex(0)];
  sigmaHigh[0] = 0.0;
  if(model) sigmaHigh[0] = model->CrossSection(couple,particle,e,cut,tmax);

  if(2 < verboseLevel) {
      G4cout << "### For material " << couple->GetMaterial()->GetName()
             << "  " << nmod
             << " models"
             << " Ecut(MeV)= " << cut/MeV
             << " Emax(MeV)= " << e/MeV
             << " nbins= "  << totBinsLambda
             << G4endl;
      G4cout << " model #0   eUp= " << e 
	     << "  sigmaUp= " << sigmaHigh[0] << G4endl;
  }


  if(nmod > 1) {

    for(j=1; j<nmod; j++) {

      e  = regModels->LowEdgeEnergy(j);
      eLow[j]   = e;
      G4int idx = regModels->ModelIndex(j); 

      sigmaLow[j] = models[idx]->CrossSection(couple,particle,e,cut,tmax);
      e  = regModels->LowEdgeEnergy(j+1);
      eLow[j+1] = e;
      sigmaHigh[j] = models[idx]->CrossSection(couple,particle,e,cut,tmax);
    }
    if(1 < verboseLevel) {
      G4cout << " model #0"  
	     << "  sigma(" << eLow[0] << ")=  " << sigmaLow[0]
	     << "  sigma(" << eLow[1] << ")=  " << sigmaHigh[0]
	     << G4endl;
    }
    for(j=1; j<nmod; j++) {
      if(sigmaLow[j] > 0.0) {
	factor[j] = (sigmaHigh[j-1]/sigmaLow[j] - 1.0)*eLow[j];
      } else  factor[j] = 0.0;
      if(1 < verboseLevel) {
        G4cout << " model #" << j 
	       << "  sigma(" << eLow[j] << ")=  " << sigmaLow[j]
	       << "  sigma(" << eLow[j+1] << ")=  " << sigmaHigh[j]
	       << "  factor= " << factor[j]/eLow[j]
	       << G4endl;
      }
    }
  }

  // Calculate lambda vector
  for(j=0; j<totBinsLambda; j++) {

    e = aVector->GetLowEdgeEnergy(j);

    // Choose a model of energy losses
    G4int k = 0;
    G4double fac = 1.0;
    if (nmod > 1 && e > eLow[1]) {
      do {
        k++;
        fac *= (1.0 + factor[k]/e);
      } while ( k+1 < nmod && e > eLow[k+1] );
    }

    model = models[regModels->ModelIndex(k)];
    G4double cross = 0.0;
    if(model) cross = model->CrossSection(couple,particle,e,cut,tmax)*fac;
    if(j==0 && startFromNull) cross = 0.0;

    if(2 < verboseLevel) {
      G4cout << "FillLambdaVector: " << j << ".   e(MeV)= " << e/MeV
	     << "  cross(1/mm)= " << cross*mm
	     << " fac= " << fac << " k= " << k 
	     << " model= " << regModels->ModelIndex(k)
	     << G4endl;
    }
    if(cross < 0.0) cross = 0.0;

    aVector->PutValue(j, cross);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::DumpModelList(G4int verb)
{
  if(verb == 0) return;
  for(G4int i=0; i<nRegions; i++) {
    G4RegionModels* r = setOfRegionModels[i];
    const G4Region* reg = r->Region();
    if(verb > 1 || nRegions > 1) {
    }
    G4int n = r->NumberOfModels();  
    if(verb > 1 || n > 0) {
      G4cout << "      ===== EM models for the G4Region  " << reg->GetName()
	     << " ======" << G4endl;;
      for(G4int j=0; j<n; j++) {
	const G4VEmModel* m = models[r->ModelIndex(j)];
	G4cout << std::setw(20);
	G4cout << m->GetName() << " :     Emin= " 
	       << std::setw(10) << G4BestUnit(r->LowEdgeEnergy(j),"Energy")
	       << "        Emax=   " 
	       << G4BestUnit(r->LowEdgeEnergy(j+1),"Energy")
	       << G4endl;
      }  
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
