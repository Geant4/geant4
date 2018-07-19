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
// $Id: G4EmModelManager.cc 106714 2017-10-20 09:38:06Z gcosmo $
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
// 08-04-08 Fixed and simplified initialisation of G4RegionModel (VI)
// 03-08-09 Create internal vectors only it is needed (VI)
// 14-07-11 Use pointer to the vector of cuts and not local copy (VI)
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
#include "G4SystemOfUnits.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4RegionModels::G4RegionModels(G4int nMod, std::vector<G4int>& indx, 
                               G4DataVector& lowE, const G4Region* reg)
{
  nModelsForRegion      = nMod;
  theListOfModelIndexes = new G4int [nModelsForRegion];
  lowKineticEnergy      = new G4double [nModelsForRegion+1];
  for (G4int i=0; i<nModelsForRegion; ++i) {
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
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4RegionStore.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4UnitsTable.hh"
#include "G4DataVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmModelManager::G4EmModelManager():
  nEmModels(0),
  nRegions(0),
  particle(0),
  verboseLevel(0)
{
  maxSubCutInRange = 0.7*mm;
  models.reserve(4);
  flucModels.reserve(4);
  regions.reserve(4);
  orderOfModels.reserve(4);
  isUsed.reserve(4);
  severalModels = true;
  fluoFlag = false;
  currRegionModel = nullptr;
  currModel  = nullptr;
  theCuts    = nullptr;
  theCutsNew = nullptr;
  theSubCuts = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmModelManager::~G4EmModelManager()
{
  verboseLevel = 0; // no verbosity at destruction
  Clear();
  delete theCutsNew; 
  delete theSubCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::Clear()
{
  if(1 < verboseLevel) {
    G4cout << "G4EmModelManager::Clear()" << G4endl;
  }
  size_t n = setOfRegionModels.size();
  if(n > 0) {
    for(size_t i=0; i<n; ++i) {
      delete setOfRegionModels[i];
      setOfRegionModels[i] = nullptr;
    }
  }
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
  isUsed.push_back(0);
  p->DefineForRegion(r);
  ++nEmModels;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::UpdateEmModel(const G4String& nam, 
                                     G4double emin, G4double emax)
{
  if (nEmModels > 0) {
    for(G4int i=0; i<nEmModels; ++i) {
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
  G4VEmModel* model = nullptr;
  if(i < nEmModels) { model = models[i]; }
  else if(verboseLevel > 0 && ver) { 
    G4cout << "G4EmModelManager::GetModel WARNING: "
           << "index " << i << " is wrong Nmodels= "
           << nEmModels;
    if(particle) { G4cout << " for " << particle->GetParticleName(); } 
    G4cout<< G4endl;
  }
  return model;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmModel* G4EmModelManager::GetRegionModel(G4int k, size_t idx)
{
  G4RegionModels* rm = setOfRegionModels[idxOfRegionModels[idx]];
  G4VEmModel* mod = models[rm->ModelIndex(k)];
  return mod;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4EmModelManager::NumberOfRegionModels(size_t idx) const
{
  G4RegionModels* rm = setOfRegionModels[idxOfRegionModels[idx]];
  return rm->NumberOfModels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4DataVector* 
G4EmModelManager::Initialise(const G4ParticleDefinition* p,
                             const G4ParticleDefinition* secondaryParticle,
                             G4double minSubRange,
                             G4int val)
{
  verboseLevel = val;
  G4String partname = p->GetParticleName();
  if(1 < verboseLevel) {
    G4cout << "G4EmModelManager::Initialise() for "
           << partname << "  Nmodels= " << nEmModels << G4endl;
  }
  // Are models defined?
  if(nEmModels < 1) {
    G4ExceptionDescription ed;
    ed << "No models found out for " << p->GetParticleName() 
       << " !";
    G4Exception("G4EmModelManager::Initialise","em0002",
                FatalException, ed);
  }

  particle = p;
  Clear(); // needed if run is not first

  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  const G4Region* world = 
    regionStore->GetRegion("DefaultRegionForTheWorld", false);

  // Identify the list of regions with different set of models
  nRegions = 1;
  std::vector<const G4Region*> setr;
  setr.push_back(world);
  G4bool isWorld = false;

  for (G4int ii=0; ii<nEmModels; ++ii) {
    const G4Region* r = regions[ii];
    if ( r == 0 || r == world) {
      isWorld = true;
      regions[ii] = world;
    } else {
      G4bool newRegion = true;
      if (nRegions>1) {
        for (G4int j=1; j<nRegions; ++j) {
          if ( r == setr[j] ) { newRegion = false; }
        }
      }
      if (newRegion) {
        setr.push_back(r);
        nRegions++;
      }
    }
  }
  // Are models defined?
  if(!isWorld) {
    G4ExceptionDescription ed;
    ed << "No models defined for the World volume for " 
       << p->GetParticleName() << " !";
    G4Exception("G4EmModelManager::Initialise","em0002",
                FatalException, ed);
  }

  G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  // prepare vectors, shortcut for the case of only 1 model
  // or only one region
  if(nRegions > 1 && nEmModels > 1) {
    idxOfRegionModels.resize(numOfCouples,0);
    setOfRegionModels.resize((size_t)nRegions,0); 
  } else {
    idxOfRegionModels.resize(1,0);
    setOfRegionModels.resize(1,0); 
  }

  std::vector<G4int>    modelAtRegion(nEmModels);
  std::vector<G4int>    modelOrd(nEmModels);
  G4DataVector          eLow(nEmModels+1);
  G4DataVector          eHigh(nEmModels);

  if(1 < verboseLevel) {
    G4cout << "    Nregions= " << nRegions 
           << "  Nmodels= " << nEmModels << G4endl;
  }

  // Order models for regions
  for (G4int reg=0; reg<nRegions; ++reg) {
    const G4Region* region = setr[reg];
    G4int n = 0;

    for (G4int ii=0; ii<nEmModels; ++ii) {

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
		 << "; tminAct= " << model->LowEnergyActivationLimit()/MeV
		 << "; tmaxAct= " << model->HighEnergyActivationLimit()/MeV
                 << G4endl;
        }
        
	static const G4double limitdelta = 0.01*eV;
        if(n > 0) {

          // extend energy range to previous models
          tmin = std::min(tmin, eHigh[n-1]);
          tmax = std::max(tmax, eLow[0]);
          //G4cout << "tmin= " << tmin << "  tmax= " 
          //           << tmax << "  ord= " << ord <<G4endl;
          // empty energy range
          if( tmax - tmin <= limitdelta) { push = false; }
          // low-energy model
          else if (tmax == eLow[0]) {
            push = false;
            insert = true;
            idx = 0;
            // resolve intersections
          } else if(tmin < eHigh[n-1]) { 
            // compare order
            for(G4int k=0; k<n; ++k) {
              // new model has higher order parameter, 
	      // so, its application area may be reduced
	      // to avoid intersections
              if(ord >= modelOrd[k]) {
                if(tmin < eHigh[k]  && tmin >= eLow[k]) { tmin = eHigh[k]; }
                if(tmax <= eHigh[k] && tmax >  eLow[k]) { tmax = eLow[k]; }
                if(tmax > eHigh[k] && tmin < eLow[k]) {
                  if(tmax - eHigh[k] > eLow[k] - tmin) { tmin = eHigh[k]; }
                  else { tmax = eLow[k]; }
                }
                if( tmax - tmin <= limitdelta) {
                  push = false;
                  break;
                }
              }
            }
	    // this model has lower order parameter than possible
	    // other models, with which there may be intersections
	    // so, appliction area of such models may be reduced
            //G4cout << "tmin= " << tmin << "  tmax= " 
            //       << tmax << "  push= " << push << " idx= " << idx <<G4endl;
	    //G4cout << "n= " << n << " eLow[0]= " <<  eLow[0] << " eLow[n-1]= " << eLow[n-1]
	    //	   << " eHigh[0]= " << eHigh[0] << " eHigh[n-1]= " << eHigh[n-1] << G4endl;

	    // insert below the first model
	    if (tmax <= eLow[0]) {
	      push = false;
	      insert = true;
	      idx = 0;
	      // resolve intersections
	    } else if(tmin < eHigh[n-1]) { 
	      // last energy interval
	      if(tmin > eLow[n-1] && tmax >= eHigh[n-1]) {
		eHigh[n-1] = tmin;
		// first energy interval
	      } else if(tmin <= eLow[0] && tmax < eHigh[0]) {
		eLow[0] = tmax;
		push = false;
		insert = true;
		idx = 0;
		// loop over all models
	      } else {
		for(G4int k=n-1; k>=0; --k) { 
		  if(tmin <= eLow[k] && tmax >= eHigh[k]) {
		    // full overlap exclude previous model
		    isUsed[modelAtRegion[k]] = 0;
		    idx = k;
		    if(k < n-1) {
		      // shift upper models and change index
		      for(G4int kk=k; kk<n-1; ++kk) {    
			modelAtRegion[kk] = modelAtRegion[kk+1];
			modelOrd[kk] = modelOrd[kk+1];
			eLow[kk] = eLow[kk+1];
			eHigh[kk] = eHigh[kk+1];
		      }
                      ++k;
		    }
		    --n;
		  } else {
		    // partially reduce previous model area
		    if(tmin <= eLow[k] && tmax > eLow[k]) { 
		      eLow[k] = tmax;
                      idx = k;
                      insert = true;
                      push = false;
		    } else if(tmin < eHigh[k] && tmax >= eHigh[k]) { 
		      eHigh[k] = tmin;
                      idx = k + 1;
                      if(idx < n) { 
			insert = true; 
			push = false;
		      }
		    } else if(tmin > eLow[k] && tmax < eHigh[k]) {
		      if(eHigh[k] - tmax > tmin - eLow[k]) {
			eLow[k] = tmax;
                        idx = k;
                        insert = true;
			push = false;
		      } else {
			eHigh[k] = tmin;
			idx = k + 1;
			if(idx < n) { 
			  insert = true; 
			  push = false;
			}
		      }
		    }
                  }
                }
              }
            }
          }
        }
	// provide space for the new model
        if(insert) {
          for(G4int k=n-1; k>=idx; --k) {    
            modelAtRegion[k+1] = modelAtRegion[k];
            modelOrd[k+1] = modelOrd[k];
            eLow[k+1]  = eLow[k];
            eHigh[k+1] = eHigh[k];
          }
        }
        //G4cout << "push= " << push << " insert= " << insert 
	//       << " idx= " << idx <<G4endl;
	// the model is added
        if (push || insert) {
          ++n;
          modelAtRegion[idx] = ii;
          modelOrd[idx] = ord;
          eLow[idx]  = tmin;
          eHigh[idx] = tmax;
          isUsed[ii] = 1;
        }
	// exclude models with zero energy range
	for(G4int k=n-1; k>=0; --k) {
	  if(eHigh[k] - eLow[k] <= limitdelta) {
	    isUsed[modelAtRegion[k]] = 0;
	    if(k < n-1) {
	      for(G4int kk=k; kk<n-1; ++kk) {    
		modelAtRegion[kk] = modelAtRegion[kk+1];
		modelOrd[kk] = modelOrd[kk+1];
		eLow[kk] = eLow[kk+1];
		eHigh[kk] = eHigh[kk+1];
	      }
	    }
	    --n;
	  }
	}
      }
    }
    eLow[0] = 0.0;
    eLow[n] = eHigh[n-1];

    if(1 < verboseLevel) {
      G4cout << "### New G4RegionModels set with " << n << " models for region <";
      if (region) { G4cout << region->GetName(); }
      G4cout << ">  Elow(MeV)= ";
      for(G4int iii=0; iii<=n; ++iii) {G4cout << eLow[iii]/MeV << " ";}
      G4cout << G4endl;
    }
    G4RegionModels* rm = new G4RegionModels(n, modelAtRegion, eLow, region);
    setOfRegionModels[reg] = rm;
    // shortcut
    if(1 == nEmModels) { break; }
  }

  currRegionModel = setOfRegionModels[0];
  currModel = models[0];

  // Access to materials and build cuts
  size_t idx = 1;
  if(secondaryParticle) {
    if( secondaryParticle == G4Gamma::Gamma() )           { idx = 0; }
    else if( secondaryParticle == G4Electron::Electron()) { idx = 1; }
    else if( secondaryParticle == G4Positron::Positron()) { idx = 2; }
    else { idx = 3; }
  }

  theCuts = 
    static_cast<const G4DataVector*>(theCoupleTable->GetEnergyCutsVector(idx));

  // for the second run the check on cuts should be repeated
  if(theCutsNew) { *theCutsNew = *theCuts; }

  if(minSubRange < 1.0) {
    if( !theSubCuts ) { theSubCuts = new G4DataVector(); }
    theSubCuts->resize(numOfCouples,DBL_MAX);
  }

  //  G4cout << "========Start define cuts" << G4endl;
  // define cut values
  for(size_t i=0; i<numOfCouples; ++i) {

    const G4MaterialCutsCouple* couple = 
      theCoupleTable->GetMaterialCutsCouple(i);
    const G4Material* material = couple->GetMaterial();
    const G4ProductionCuts* pcuts = couple->GetProductionCuts();
 
    G4int reg = 0;
    if(nRegions > 1 && nEmModels > 1) {
      reg = nRegions;
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
      do {--reg;} while (reg>0 && pcuts != (setr[reg]->GetProductionCuts()));
      idxOfRegionModels[i] = reg;
    }
    if(1 < verboseLevel) {
      G4cout << "G4EmModelManager::Initialise() for "
             << material->GetName()
             << " indexOfCouple= " << i
             << " indexOfRegion= " << reg
             << G4endl;
    }

    G4double cut = (*theCuts)[i]; 
    if(secondaryParticle) {

      // compute subcut 
      if( cut < DBL_MAX && minSubRange < 1.0) {
        G4double subcut = minSubRange*cut;
        G4double rcut = std::min(minSubRange*pcuts->GetProductionCut(idx), 
                                 maxSubCutInRange);
        G4double tcutmax = 
          theCoupleTable->ConvertRangeToEnergy(secondaryParticle,
                                               material,rcut);
        if(tcutmax < subcut) { subcut = tcutmax; }
        (*theSubCuts)[i] = subcut;
      }

      // note that idxOfRegionModels[] not always filled
      G4int inn = 0;
      G4int nnm = 1;
      if(nRegions > 1 && nEmModels > 1) { 
        inn = idxOfRegionModels[i];
      }
      // check cuts and introduce upper limits
      //G4cout << "idx= " << i << " cut(keV)= " << cut/keV << G4endl;
      currRegionModel = setOfRegionModels[inn];
      nnm = currRegionModel->NumberOfModels();
     
      //G4cout << "idx= " << i << " Nmod= " << nnm << G4endl; 

      for(G4int jj=0; jj<nnm; ++jj) {
        //G4cout << "jj= " << jj << "  modidx= " 
        //       << currRegionModel->ModelIndex(jj) << G4endl;
        currModel = models[currRegionModel->ModelIndex(jj)];
        G4double cutlim = currModel->MinEnergyCut(particle,couple);
        if(cutlim > cut) {
          if(!theCutsNew) { theCutsNew = new G4DataVector(*theCuts); }
          (*theCutsNew)[i] = cutlim;
          /*
          G4cout << "### " << partname << " energy loss model in " 
                 << material->GetName() 
                 << "  Cut was changed from " << cut/keV << " keV to "
                 << cutlim/keV << " keV " << " due to " 
                 << currModel->GetName() << G4endl;
          */
        }
      } 
    }
  }
  if(theCutsNew) { theCuts = theCutsNew; }

  // initialize models
  G4int nn = 0;
  severalModels = true;
  for(G4int jj=0; jj<nEmModels; ++jj) {
    if(1 == isUsed[jj]) {
      ++nn;
      currModel = models[jj];
      currModel->Initialise(particle, *theCuts);
      if(flucModels[jj]) { flucModels[jj]->InitialiseMe(particle); }
    }
  }
  if(1 == nn) { severalModels = false; }

  if(1 < verboseLevel) {
    G4cout << "G4EmModelManager for " << partname 
           << " is initialised; nRegions=  " << nRegions
           << " severalModels: " << severalModels 
           << G4endl;
  }

  return theCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::FillDEDXVector(G4PhysicsVector* aVector,
                                      const G4MaterialCutsCouple* couple,
                                      G4EmTableType tType)
{
  size_t i = couple->GetIndex();
  G4double cut  = (*theCuts)[i];
  G4double emin = 0.0;

  if(fTotal == tType) { cut = DBL_MAX; }
  else if(fSubRestricted == tType) {
    emin = cut;
    if(theSubCuts) { emin = (*theSubCuts)[i]; }
  }

  if(1 < verboseLevel) {
    G4cout << "G4EmModelManager::FillDEDXVector() for "
           << couple->GetMaterial()->GetName()
           << "  cut(MeV)= " << cut
           << "  emin(MeV)= " << emin
           << "  Type " << tType
           << "  for " << particle->GetParticleName()
           << G4endl;
  }

  G4int reg  = 0;
  if(nRegions > 1 && nEmModels > 1) { reg = idxOfRegionModels[i]; }
  const G4RegionModels* regModels = setOfRegionModels[reg];
  G4int nmod = regModels->NumberOfModels();

  // Calculate energy losses vector

  //G4cout << "nmod= " << nmod << G4endl;
  size_t totBinsLoss = aVector->GetVectorLength();
  G4double del = 0.0;
  G4int    k0  = 0;

  for(size_t j=0; j<totBinsLoss; ++j) {

    G4double e = aVector->Energy(j);

    // Choose a model of energy losses
    G4int k = 0;
    if (nmod > 1) {
      k = nmod;
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
      do {--k;} while (k>0 && e <= regModels->LowEdgeEnergy(k));
      //G4cout << "k= " << k << G4endl;
      if(k > 0 && k != k0) {
        k0 = k;
        G4double elow = regModels->LowEdgeEnergy(k);
        G4double dedx1 = ComputeDEDX(models[regModels->ModelIndex(k-1)],
                                     couple,elow,cut,emin);
        G4double dedx2 = ComputeDEDX(models[regModels->ModelIndex(k)],
                                     couple,elow,cut,emin);
        del = 0.0;
        if(dedx2 > 0.0) { del = (dedx1/dedx2 - 1.0)*elow; }
        //G4cout << "elow= " << elow 
        //       << " dedx1= " << dedx1 << " dedx2= " << dedx2 << G4endl;
      }
    }
    G4double dedx = 
      ComputeDEDX(models[regModels->ModelIndex(k)],couple,e,cut,emin);
    dedx *= (1.0 + del/e); 

    if(2 < verboseLevel) {
      G4cout << "Material= " << couple->GetMaterial()->GetName()
             << "   E(MeV)= " << e/MeV
             << "  dEdx(MeV/mm)= " << dedx*mm/MeV
             << "  del= " << del*mm/MeV<< " k= " << k 
             << " modelIdx= " << regModels->ModelIndex(k)
             << G4endl;
    }
    if(dedx < 0.0) { dedx = 0.0; }
    aVector->PutValue(j, dedx);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::FillLambdaVector(G4PhysicsVector* aVector,
                                        const G4MaterialCutsCouple* couple,
                                        G4bool startFromNull,
                                        G4EmTableType tType)
{
  size_t i = couple->GetIndex();
  G4double cut  = (*theCuts)[i];
  G4double tmax = DBL_MAX;
  if (fSubRestricted == tType) {
    tmax = cut;
    if(theSubCuts) { cut  = (*theSubCuts)[i]; }
  }

  G4int reg  = 0;
  if(nRegions > 1 && nEmModels > 1) { reg  = idxOfRegionModels[i]; }
  const G4RegionModels* regModels = setOfRegionModels[reg];
  G4int nmod = regModels->NumberOfModels();
  if(1 < verboseLevel) {
    G4cout << "G4EmModelManager::FillLambdaVector() for "
           << particle->GetParticleName()
           << " in " << couple->GetMaterial()->GetName()
           << " Emin(MeV)= " << aVector->Energy(0)
           << " Emax(MeV)= " << aVector->GetMaxEnergy()
           << " cut= " << cut
           << " Type " << tType   
           << " nmod= " << nmod
           << " theSubCuts " << theSubCuts
           << G4endl;
  }

  // Calculate lambda vector
  size_t totBinsLambda = aVector->GetVectorLength();
  G4double del = 0.0;
  G4int    k0  = 0;
  G4int     k  = 0;
  G4VEmModel* mod = models[regModels->ModelIndex(0)]; 
  for(size_t j=0; j<totBinsLambda; ++j) {

    G4double e = aVector->Energy(j);

    // Choose a model 
    if (nmod > 1) {
      k = nmod;
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
      do {--k;} while (k>0 && e <= regModels->LowEdgeEnergy(k));
      if(k > 0 && k != k0) {
        k0 = k;
        G4double elow = regModels->LowEdgeEnergy(k);
        G4VEmModel* mod1 = models[regModels->ModelIndex(k-1)]; 
        G4double xs1  = mod1->CrossSection(couple,particle,elow,cut,tmax);
        mod = models[regModels->ModelIndex(k)]; 
        G4double xs2 = mod->CrossSection(couple,particle,elow,cut,tmax);
        del = 0.0;
        if(xs2 > 0.0) { del = (xs1/xs2 - 1.0)*elow; }
        //G4cout << "New model k=" << k << " E(MeV)= " << e/MeV 
        //       << " Elow(MeV)= " << elow/MeV << " del= " << del << G4endl;
      }
    }
    G4double cross = mod->CrossSection(couple,particle,e,cut,tmax);
    cross *= (1.0 + del/e); 
    if(fIsCrossSectionPrim == tType) { cross *= e; }
    
    if(j==0 && startFromNull) { cross = 0.0; }

    if(2 < verboseLevel) {
      G4cout << "FillLambdaVector: " << j << ".   e(MeV)= " << e/MeV
             << "  cross(1/mm)= " << cross*mm
             << " del= " << del*mm << " k= " << k 
             << " modelIdx= " << regModels->ModelIndex(k)
             << G4endl;
    }
    if(cross < 0.0) { cross = 0.0; }
    aVector->PutValue(j, cross);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::DumpModelList(std::ostream& out, G4int verb)
{
  DumpModelList(out, verb, G4String("\n"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EmModelManager::DumpModelList(std::ostream& out, G4int verb,
                                     G4String endOfLine)
{
  if(verb == 0) { return; }
  for(G4int i=0; i<nRegions; ++i) {
    G4RegionModels* r = setOfRegionModels[i];
    const G4Region* reg = r->Region();
    G4int n = r->NumberOfModels();  
    if(n > 0) {
      out << "      ===== EM models for the G4Region  " << reg->GetName()
	  << " ======" << endOfLine;
      for(G4int j=0; j<n; ++j) {
        G4VEmModel* model = models[r->ModelIndex(j)];
        G4double emin = 
          std::max(r->LowEdgeEnergy(j),model->LowEnergyActivationLimit());
        G4double emax = 
          std::min(r->LowEdgeEnergy(j+1),model->HighEnergyActivationLimit());
        if(emax > emin) {
	  out << std::setw(20);
	  out << model->GetName() << " :  Emin= " 
	      << std::setw(8) << G4BestUnit(emin,"Energy")
	      << "   Emax= " 
	      << std::setw(8) << G4BestUnit(emax,"Energy");
	  G4PhysicsTable* table = model->GetCrossSectionTable();
	  if(table) {
	    size_t kk = table->size();
	    for(size_t k=0; k<kk; ++k) {
	      G4PhysicsVector* v = (*table)[k];
	      if(v) {
		G4int nn = v->GetVectorLength() - 1;
		out << "  Table with " << nn << " bins Emin= "
		    << std::setw(6) << G4BestUnit(v->Energy(0),"Energy")
		    << "   Emax= " 
		    << std::setw(6) << G4BestUnit(v->Energy(nn),"Energy");
		break;
	      }
	    }
          }
	  G4VEmAngularDistribution* an = model->GetAngularDistribution();
	  if(an) { out << "   " << an->GetName(); }
	  if(fluoFlag && model->DeexcitationFlag()) { 
	    out << "  FluoActive"; 
	  }
	  out << endOfLine;
	}
      }  
    }
    if(1 == nEmModels) { break; }
  }
  if(theCutsNew) {
    out << "      ===== Limit on energy threshold has been applied " 
	<< endOfLine;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
