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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4LossTableBuilder
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 23-01-03 V.Ivanchenko Cut per region
// 21-07-04 V.Ivanchenko Fix problem of range for dedx=0
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivanchenko)
// 07-12-04 Fix of BuildDEDX table (V.Ivanchenko)
// 27-03-06 Add bool options isIonisation (V.Ivanchenko)
// 16-01-07 Fill new (not old) DEDX table (V.Ivanchenko)
// 12-02-07 Use G4LPhysicsFreeVector for the inverse range table (V.Ivanchenko)
// 24-06-09 Removed hidden bin in G4PhysicsVector (V.Ivanchenko)
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LossTableBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4VEmModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"

std::vector<G4double>* G4LossTableBuilder::theDensityFactor = nullptr;
std::vector<G4int>*    G4LossTableBuilder::theDensityIdx = nullptr;
std::vector<G4bool>*   G4LossTableBuilder::theFlag = nullptr;
#ifdef G4MULTITHREADED
  G4Mutex G4LossTableBuilder::ltbMutex = G4MUTEX_INITIALIZER;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LossTableBuilder::G4LossTableBuilder(G4bool master) : isMaster(master) 
{
  theParameters = G4EmParameters::Instance();
  if(nullptr == theFlag) {
#ifdef G4MULTITHREADED
    G4MUTEXLOCK(&ltbMutex);
    if(nullptr == theFlag) {
#endif
      if(!isMaster) {
	G4ExceptionDescription ed;
	ed << "Initialisation called from a worker thread ";
	G4Exception("G4LossTableBuilder: ", "em0001",
		    JustWarning, ed);
      }
      theDensityFactor = new std::vector<G4double>;
      theDensityIdx = new std::vector<G4int>;
      theFlag = new std::vector<G4bool>;
#ifdef G4MULTITHREADED
    }
    G4MUTEXUNLOCK(&ltbMutex);
#endif
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LossTableBuilder::~G4LossTableBuilder() 
{
  if(isMaster) {
    delete theDensityFactor;
    delete theDensityIdx;
    delete theFlag;
    theDensityFactor = nullptr;
    theDensityIdx = nullptr;
    theFlag = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::vector<G4int>* G4LossTableBuilder::GetCoupleIndexes() const
{
  return theDensityIdx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::vector<G4double>* G4LossTableBuilder::GetDensityFactors() const
{
  return theDensityFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4LossTableBuilder::GetFlag(std::size_t idx)
{
  if(theFlag->empty()) { InitialiseBaseMaterials(); }
  return (idx < theFlag->size()) ? (*theFlag)[idx] : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4LossTableBuilder::GetBaseMaterialFlag()
{
  if(theFlag->empty()) { InitialiseBaseMaterials(); }
  return baseMatFlag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LossTableBuilder::BuildDEDXTable(G4PhysicsTable* dedxTable,
                                   const std::vector<G4PhysicsTable*>& list)
{
  InitialiseBaseMaterials(dedxTable);
  std::size_t n_processes = list.size();
  if(1 >= n_processes) { return; }

  std::size_t nCouples = dedxTable->size();
  //G4cout << "Nproc= " << n_processes << " nCouples=" << nCouples << " Nv= " 
  //	 << dedxTable->size() << G4endl;
  if(0 >= nCouples) { return; }

  for (std::size_t i=0; i<nCouples; ++i) {
    auto pv0 = static_cast<G4PhysicsLogVector*>((*(list[0]))[i]);
    if(pv0 == nullptr) { continue; } 
    std::size_t npoints = pv0->GetVectorLength();
    auto pv = new G4PhysicsLogVector(*pv0);
    for (std::size_t j=0; j<npoints; ++j) {
      G4double dedx = 0.0;
      for (std::size_t k=0; k<n_processes; ++k) {
	const G4PhysicsVector* pv1 = (*(list[k]))[i];
	dedx += (*pv1)[j];
      }
      pv->PutValue(j, dedx);
    }
    if(splineFlag) { pv->FillSecondDerivatives(); }
    G4PhysicsTableHelper::SetPhysicsVector(dedxTable, i, pv);
  }
  //G4cout << "### G4LossTableBuilder::BuildDEDXTable " << G4endl; 
  //G4cout << *dedxTable << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::BuildRangeTable(const G4PhysicsTable* dedxTable,
                                         G4PhysicsTable* rangeTable)
// Build range table from the energy loss table
{
  //G4cout << "### G4LossTableBuilder::BuildRangeTable: DEDX table" << G4endl; 
  //G4cout << *const_cast<G4PhysicsTable*>(dedxTable) << G4endl;
  const std::size_t nCouples = dedxTable->size();
  if(0 >= nCouples) { return; }

  const std::size_t n = 100;
  const G4double del = 1.0/(G4double)n;

  for (std::size_t i=0; i<nCouples; ++i) {
    auto pv = static_cast<G4PhysicsLogVector*>((*dedxTable)[i]);
    if((pv == nullptr) || (isBaseMatActive && !(*theFlag)[i])) { continue; } 
    std::size_t npoints = pv->GetVectorLength();
    std::size_t bin0    = 0;
    G4double elow  = pv->Energy(0);
    G4double ehigh = pv->Energy(npoints-1);
    G4double dedx1 = (*pv)[0];

    //G4cout << "i= " << i << "npoints= " << npoints << " dedx1= " 
    //<< dedx1 << G4endl;

    // protection for specific cases dedx=0
    if(dedx1 == 0.0) {
      for (std::size_t k=1; k<npoints; ++k) {
	++bin0;
	elow  = pv->Energy(k);
	dedx1 = (*pv)[k];
	if(dedx1 > 0.0) { break; }
      }
      npoints -= bin0;
    }
    //G4cout<<"New Range vector" << G4endl;
    //G4cout<<"nbins= "<<npoints-1<<" elow= "<<elow<<" ehigh= "<<ehigh
    //            <<" bin0= " << bin0 <<G4endl;

    // initialisation of a new vector
    if(npoints < 3) { npoints = 3; }

    delete (*rangeTable)[i];
    G4PhysicsLogVector* v;
    if(0 == bin0) { v = new G4PhysicsLogVector(*pv); }
    else { v = new G4PhysicsLogVector(elow, ehigh, npoints-1, splineFlag); }

    // assumed dedx proportional to beta
    G4double energy1 = v->Energy(0);
    G4double range   = 2.*energy1/dedx1;
    //G4cout << "range0= " << range << G4endl;
    v->PutValue(0,range);

    for (std::size_t j=1; j<npoints; ++j) {

      G4double energy2 = v->Energy(j);
      G4double de      = (energy2 - energy1) * del;
      G4double energy  = energy2 + de*0.5;
      G4double sum = 0.0;
      std::size_t idx = j - 1;
      //G4cout << "j= " << j << " e1= " << energy1 << " e2= " << energy2 
      //       << " n= " << n << G4endl;
      for (std::size_t k=0; k<n; ++k) {
	energy -= de;
	dedx1 = pv->Value(energy, idx);
	if(dedx1 > 0.0) { sum += de/dedx1; }
      }
      range += sum;
      v->PutValue(j,range);
      energy1 = energy2;
    }
    if(splineFlag) { v->FillSecondDerivatives(); }
    G4PhysicsTableHelper::SetPhysicsVector(rangeTable, i, v);
  }
  //G4cout << "### Range table" << G4endl; 
  //G4cout << *rangeTable << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LossTableBuilder::BuildInverseRangeTable(const G4PhysicsTable* rangeTable,
                                           G4PhysicsTable* invRangeTable)
// Build inverse range table from the energy loss table
{
  std::size_t nCouples = rangeTable->size();
  if(0 >= nCouples) { return; }

  for (std::size_t i=0; i<nCouples; ++i) {
    G4PhysicsVector* pv = (*rangeTable)[i];
    if((pv == nullptr) || (isBaseMatActive && !(*theFlag)[i])) { continue; } 
    std::size_t npoints = pv->GetVectorLength();
      
    delete (*invRangeTable)[i];
    auto v = new G4PhysicsFreeVector(npoints,splineFlag);

    for (std::size_t j=0; j<npoints; ++j) {
      G4double e  = pv->Energy(j);
      G4double r  = (*pv)[j];
      v->PutValues(j,r,e);
    }
    if(splineFlag) { v->FillSecondDerivatives(); }

    G4PhysicsTableHelper::SetPhysicsVector(invRangeTable, i, v);
  }
  //G4cout << "### Inverse range table" << G4endl; 
  //G4cout << *invRangeTable << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::InitialiseBaseMaterials(const G4PhysicsTable* table)
{
  if(!isMaster) { return; }
  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t nCouples = theCoupleTable->GetTableSize();
  std::size_t nFlags = theFlag->size();
  /*
  G4cout << "### InitialiseBaseMaterials: nCouples=" << nCouples
	 << " nFlags=" << nFlags << " isInit:" << isInitialized
	 << " baseMat:" << baseMatFlag << G4endl; 
  */
  // define base material flag
  if(isBaseMatActive && !baseMatFlag) {
    for(G4int i=0; i<(G4int)nCouples; ++i) {
      if(nullptr != theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial()->GetBaseMaterial()) {
	baseMatFlag = true;
        isInitialized = false;
	break;
      }
    }
  }

  if(nFlags != nCouples) { isInitialized = false; }
  if(isInitialized) { return; }

  // reserve memory
  theFlag->resize(nCouples, true);
  if(nullptr == table) { return; }

  if(baseMatFlag) {
    theDensityFactor->resize(nCouples,1.0);
    theDensityIdx->resize(nCouples);
  }

  // define default flag and index of used material cut couple
  for(G4int i=0; i<(G4int)nCouples; ++i) {
    (*theFlag)[i] = table->GetFlag(i); 
    if(baseMatFlag) { (*theDensityIdx)[i] = i; }
  }
  isInitialized = true;
  if(baseMatFlag) {
    // use base materials
    for(G4int i=0; i<(G4int)nCouples; ++i) { 
      // base material is needed only for a couple which is not
      // initialised and for which tables will be computed
      auto couple = theCoupleTable->GetMaterialCutsCouple(i);
      auto pcuts = couple->GetProductionCuts();
      auto mat  = couple->GetMaterial();
      auto bmat = mat->GetBaseMaterial();

      // base material exists - find it and check if it can be reused
      if(nullptr != bmat) {
	for(G4int j=0; j<(G4int)nCouples; ++j) {
	  if(j == i) { continue; }
	  auto bcouple = theCoupleTable->GetMaterialCutsCouple(j);

	  if(bcouple->GetMaterial() == bmat && 
	     bcouple->GetProductionCuts() == pcuts) {

	    // based couple exist in the same region
	    (*theDensityFactor)[i] = mat->GetDensity()/bmat->GetDensity();
	    (*theDensityIdx)[i] = j;
	    (*theFlag)[i] = false;

	    // ensure that there will no double initialisation
	    (*theDensityFactor)[j] = 1.0;
	    (*theDensityIdx)[j] = j;
	    (*theFlag)[j] = true;
	    break;
	  } 
	}
      }
    }
  }
  /*
  G4cout << "### G4LossTableBuilder::InitialiseBaseMaterials: flag=" 
	 << baseMatFlag << G4endl;
  for(std::size_t i=0; i<nCouples; ++i) {
    G4cout << "CoupleIdx=" << i << "  Flag= " << (*theFlag)[i] << "  "
           << theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial()->GetName()
	   << "  TableFlag= " << table->GetFlag(i)
           << "  " << (*table)[i]
	   << G4endl;
  }
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* 
G4LossTableBuilder::BuildTableForModel(G4PhysicsTable* aTable, 
                                       G4VEmModel* model, 
                                       const G4ParticleDefinition* part,
                                       G4double emin, G4double emax,
                                       G4bool spline)
{
  // check input
  G4PhysicsTable* table = G4PhysicsTableHelper::PreparePhysicsTable(aTable);
  if(nullptr == table) { return table; }
  if(emin >= emax) { 
    table->clearAndDestroy();
    delete table;
    table = nullptr;
    return table; 
  }
  InitialiseBaseMaterials(table);
  G4int nbins = theParameters->NumberOfBinsPerDecade();

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  std::size_t numOfCouples = theCoupleTable->GetTableSize();

  G4PhysicsLogVector* aVector = nullptr;

  for(G4int i=0; i<(G4int)numOfCouples; ++i) {
    if ((*theFlag)[i]) {

      // create physics vector and fill it
      auto couple = theCoupleTable->GetMaterialCutsCouple(i);
      delete (*table)[i];

      // if start from zero then change the scale

      const G4Material* mat = couple->GetMaterial();

      G4double tmin = std::max(emin,model->MinPrimaryEnergy(mat,part));
      if(0.0 >= tmin) { tmin = CLHEP::eV; }
      G4int n = nbins;

      if(tmin >= emax) {
        aVector = nullptr;
      } else {
        n *= G4lrint(std::log10(emax/tmin));
        n = std::max(n, 3);
        aVector = new G4PhysicsLogVector(tmin, emax, n, spline);
      }

      if(nullptr != aVector) {
        //G4cout << part->GetParticleName() << " in " << mat->GetName() 
	//       << " tmin= " << tmin << G4endl;
        for(G4int j=0; j<=n; ++j) {
          aVector->PutValue(j, model->Value(couple, part, 
                                            aVector->Energy(j)));
        }
        if(spline) { aVector->FillSecondDerivatives(); }
      }
      G4PhysicsTableHelper::SetPhysicsVector(table, i, aVector);
    }
  }
  /*
  G4cout << "G4LossTableBuilder::BuildTableForModel done for "
         << part->GetParticleName() << " and "<< model->GetName()
         << "  " << table << G4endl;
  */
  //G4cout << *table << G4endl;
  return table; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
