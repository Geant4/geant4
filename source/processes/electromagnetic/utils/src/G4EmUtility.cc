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
// Geant4 class G4EmUtility
//
// Author V.Ivanchenko 14.03.2022
//

#include "G4EmUtility.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCutsTable.hh"
#include "G4VEmProcess.hh"
#include "G4EmParameters.hh"
#include "G4PhysicsVector.hh"
#include "Randomize.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static const G4double g4log10 = G4Log(10.);

const G4Region* 
G4EmUtility::FindRegion(const G4String& regionName, const G4int verbose)
{
  const G4Region* reg = nullptr;
  G4RegionStore* regStore = G4RegionStore::GetInstance();
  G4String r = regionName;
  if(r == "") { r = "DefaultRegionForTheWorld"; }
  reg = regStore->GetRegion(r, true); 
  if(nullptr == reg && verbose > 0) {
    G4cout << "### G4EmUtility WARNING: fails to find a region <"
           << r << G4endl;
  } else if(verbose > 1) {
    G4cout << "### G4EmUtility finds out G4Region <" << r << ">" 
           << G4endl;
  } 
  return reg;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const G4Element* G4EmUtility::SampleRandomElement(const G4Material* mat)
{
  const G4Element* elm = mat->GetElement(0);
  std::size_t nElements = mat->GetNumberOfElements();
  if(1 < nElements) {
    G4double x = mat->GetTotNbOfElectPerVolume()*G4UniformRand();
    const G4double* y = mat->GetVecNbOfAtomsPerVolume();
    for(std::size_t i=0; i<nElements; ++i) {
      elm = mat->GetElement((G4int)i);
      x -= y[i]*elm->GetZ();
      if(x <= 0.0) { break; }
    }
  }
  return elm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const G4Isotope* G4EmUtility::SampleRandomIsotope(const G4Element* elm)
{
  const std::size_t ni = elm->GetNumberOfIsotopes();
  const G4Isotope* iso = elm->GetIsotope(0);
  if(ni > 1) {
    const G4double* ab = elm->GetRelativeAbundanceVector();
    G4double x = G4UniformRand();
    for(std::size_t idx=0; idx<ni; ++idx) {
      x -= ab[idx];
      if (x <= 0.0) { 
	iso = elm->GetIsotope((G4int)idx);
	break; 
      }
    }
  }
  return iso;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::vector<G4double>* G4EmUtility::FindCrossSectionMax(G4PhysicsTable* p)
{
  std::vector<G4double>* ptr = nullptr;
  if(nullptr == p) { return ptr; }

  const std::size_t n = p->length();
  ptr = new std::vector<G4double>;
  ptr->resize(n, DBL_MAX);

  G4bool isPeak = false;
  G4double e, ss, ee, xs;

  // first loop on existing vectors
  for (std::size_t i=0; i<n; ++i) {
    const G4PhysicsVector* pv = (*p)[i];
    xs = ee = 0.0;
    if(nullptr != pv) {
      G4int nb = (G4int)pv->GetVectorLength();
      for (G4int j=0; j<nb; ++j) {
	e = pv->Energy(j);
	ss = (*pv)(j);
	if(ss >= xs) {
	  xs = ss;
	  ee = e;
	  continue;
	} else {
	  isPeak = true;
	  (*ptr)[i] = ee;
	  break;
	}
      }
    }
  }

  // there is no peak for any material
  if(!isPeak) {
    delete ptr;
    ptr = nullptr;
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::vector<G4double>* 
G4EmUtility::FindCrossSectionMax(G4VDiscreteProcess* p,
                                 const G4ParticleDefinition* part)
{
  std::vector<G4double>* ptr = nullptr;
  if(nullptr == p || nullptr == part) { return ptr; }
  /*
  G4cout << "G4EmUtility::FindCrossSectionMax for "
	 << p->GetProcessName() << " and " << part->GetParticleName() << G4endl;
  */
  G4EmParameters* theParameters = G4EmParameters::Instance();
  G4double tmin = theParameters->MinKinEnergy();
  G4double tmax = theParameters->MaxKinEnergy();

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t n = theCoupleTable->GetTableSize();
  ptr = new std::vector<G4double>;
  ptr->resize(n, DBL_MAX);

  G4bool isPeak = false;
  G4double scale = theParameters->NumberOfBinsPerDecade()/g4log10;

  G4double e, sig, ee, x, sm, em, emin, emax;

  // first loop on existing vectors
  for (std::size_t i=0; i<n; ++i) {
    auto couple = theCoupleTable->GetMaterialCutsCouple((G4int)i);
    emin = std::max(p->MinPrimaryEnergy(part, couple->GetMaterial()), tmin);
    emax = std::max(tmax, 2*emin);
    ee = G4Log(emax/emin);

    G4int nbin = G4lrint(ee*scale);
    if(nbin < 4) { nbin = 4; }
    x = G4Exp(ee/nbin);
    sm = 0.0;
    em = 0.0;
    e = emin;
    for(G4int j=0; j<=nbin; ++j) {
      sig = p->GetCrossSection(e, couple);
      if(sig >= sm) {
	em = e;
	sm = sig;
	e = (j+1 < nbin) ? e*x : emax;
      } else {
	isPeak = true;
	(*ptr)[i] = em;
	break;
      }
    }
    //G4cout << i << ".  em=" << em << " sm=" << sm << G4endl;
  }
  // there is no peak for any couple
  if(!isPeak) {
    delete ptr;
    ptr = nullptr;
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

std::vector<G4TwoPeaksXS*>* 
G4EmUtility::FillPeaksStructure(G4PhysicsTable* p, G4LossTableBuilder* bld)
{
  std::vector<G4TwoPeaksXS*>* ptr = nullptr;
  if(nullptr == p) { return ptr; }

  const G4int n = (G4int)p->length();
  ptr = new std::vector<G4TwoPeaksXS*>;
  ptr->resize(n, nullptr);

  G4double e, ss, xs, ee;
  G4double e1peak, e1deep, e2peak, e2deep, e3peak;
  G4bool isDeep = false;

  // first loop on existing vectors
  for (G4int i=0; i<n; ++i) {
    const G4PhysicsVector* pv = (*p)[i];
    ee = xs = 0.0;
    e1peak = e1deep = e2peak = e2deep = e3peak = DBL_MAX;
    if(nullptr != pv) {
      G4int nb = (G4int)pv->GetVectorLength();
      for (G4int j=0; j<nb; ++j) {
	e = pv->Energy(j);
	ss = (*pv)(j);
	// find out 1st peak
	if(e1peak == DBL_MAX) {
	  if(ss >= xs) {
	    xs = ss;
	    ee = e;
	    continue;
	  } else {
	    e1peak = ee;
	  }
	}
	// find out the deep
	if(e1deep == DBL_MAX) {
	  if(ss <= xs) {
	    xs = ss;
	    ee = e;
	    continue;
	  } else {
	    e1deep = ee;
	    isDeep = true;
	  }
	}
	// find out 2nd peak
	if(e2peak == DBL_MAX) {
	  if(ss >= xs) {
	    xs = ss;
	    ee = e;
	    continue;
	  } else {
	    e2peak = ee;
	  }
	}
	if(e2deep == DBL_MAX) {
	  if(ss <= xs) {
	    xs = ss;
	    ee = e;
	    continue;
	  } else {
	    e2deep = ee;
	    break;
	  }
	}
	// find out 3d peak
	if(e3peak == DBL_MAX) {
	  if(ss >= xs) {
	    xs = ss;
	    ee = e;
	    continue;
	  } else {
	    e3peak = ee;
	  }
	}
      }
    }
    G4TwoPeaksXS* x = (*ptr)[i];
    if(nullptr == x) { 
      x = new G4TwoPeaksXS(); 
      (*ptr)[i] = x;
    }
    x->e1peak = e1peak;
    x->e1deep = e1deep;
    x->e2peak = e2peak;
    x->e2deep = e2deep;
    x->e3peak = e3peak;
  }
  // case of no 1st peak in all vectors
  if(!isDeep) {
    delete ptr;
    ptr = nullptr;
    return ptr;
  }
  // check base particles
  if(!bld->GetBaseMaterialFlag()) { return ptr; }

  auto theDensityIdx = bld->GetCoupleIndexes();
  // second loop using base materials
  for (G4int i=0; i<n; ++i) {
    const G4PhysicsVector* pv = (*p)[i];
    if (nullptr == pv) {
      G4int j = (*theDensityIdx)[i];
      if(j == i) { continue; }
      G4TwoPeaksXS* x = (*ptr)[i];
      G4TwoPeaksXS* y = (*ptr)[j];
      if(nullptr == x) { 
	x = new G4TwoPeaksXS(); 
	(*ptr)[i] = x;
      }
      x->e1peak = y->e1peak;
      x->e1deep = y->e1deep;
      x->e2peak = y->e2peak;
      x->e2deep = y->e2deep;
      x->e3peak = y->e3peak;
    }
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4EmUtility::InitialiseElementSelectors(G4VEmModel* mod,
					     const G4ParticleDefinition* part,
					     const G4DataVector& cuts,
                                             const G4double elow,
                                             const G4double ehigh)
{
  // using spline for element selectors should be investigated in details
  // because small number of points may provide biased results
  // large number of points requires significant increase of memory
  G4bool spline = false;

  G4int nbinsPerDec = G4EmParameters::Instance()->NumberOfBinsPerDecade();

  G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numOfCouples = theCoupleTable->GetTableSize();

  // prepare vector
  auto elmSelectors = mod->GetElementSelectors();
  if(nullptr == elmSelectors) {
    elmSelectors = new std::vector<G4EmElementSelector*>;
  }
  std::size_t nSelectors = elmSelectors->size();
  if(numOfCouples > nSelectors) { 
    for(std::size_t i=nSelectors; i<numOfCouples; ++i) { 
      elmSelectors->push_back(nullptr); 
    }
    nSelectors = numOfCouples;
  }

  // initialise vector
  for(std::size_t i=0; i<numOfCouples; ++i) {

    // no need in element selectors for infinite cuts
    if(cuts[i] == DBL_MAX) { continue; }
   
    auto couple = theCoupleTable->GetMaterialCutsCouple((G4int)i); 
    auto mat = couple->GetMaterial();
    mod->SetCurrentCouple(couple);

    // selector already exist then delete
    delete (*elmSelectors)[i];

    G4double emin = std::max(elow, mod->MinPrimaryEnergy(mat, part, cuts[i]));
    G4double emax = std::max(ehigh, 10*emin);
    static const G4double invlog106 = 1.0/(6*G4Log(10.));
    G4int nbins = G4lrint(nbinsPerDec*G4Log(emax/emin)*invlog106);
    nbins = std::max(nbins, 3);

    (*elmSelectors)[i] = new G4EmElementSelector(mod,mat,nbins,
						 emin,emax,spline);
    ((*elmSelectors)[i])->Initialise(part, cuts[i]);
    /*      
      G4cout << "G4VEmModel::InitialiseElmSelectors i= " << i 
             << "  "  << part->GetParticleName() 
             << " for " << mod->GetName() << "  cut= " << cuts[i] 
             << "  " << (*elmSelectors)[i] << G4endl;      
      ((*elmSelectors)[i])->Dump(part);
    */
  }
  mod->SetElementSelectors(elmSelectors);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
