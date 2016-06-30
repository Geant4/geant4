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
// $Id: G4LossTableBuilder.cc 96088 2016-03-14 16:03:38Z gcosmo $
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
#include "G4LPhysicsFreeVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4VEmModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LossTableBuilder::G4LossTableBuilder() 
{
  theParameters = G4EmParameters::Instance();
  splineFlag = true;
  isInitialized = false;

  theDensityFactor = new std::vector<G4double>;
  theDensityIdx = new std::vector<G4int>;
  theFlag = new std::vector<G4bool>;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LossTableBuilder::~G4LossTableBuilder() 
{
  delete theDensityFactor;
  delete theDensityIdx;
  delete theFlag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LossTableBuilder::BuildDEDXTable(G4PhysicsTable* dedxTable,
                                   const std::vector<G4PhysicsTable*>& list)
{
  size_t n_processes = list.size();
  //G4cout << "Nproc= " << n_processes << " Ncoup= " 
  //<< dedxTable->size() << G4endl;
  if(1 >= n_processes) { return; }

  size_t nCouples = dedxTable->size();
  if(0 >= nCouples) { return; }

  for (size_t i=0; i<nCouples; ++i) {
    //    if ((*theFlag)[i]) {
    G4PhysicsLogVector* pv0 = 
      static_cast<G4PhysicsLogVector*>((*(list[0]))[i]);
    if(pv0) {
      size_t npoints = pv0->GetVectorLength();
      G4PhysicsLogVector* pv = new G4PhysicsLogVector(*pv0);
      pv->SetSpline(splineFlag);
      for (size_t j=0; j<npoints; ++j) {
        G4double dedx = 0.0;
        for (size_t k=0; k<n_processes; ++k) {
          G4PhysicsVector* pv1   = (*(list[k]))[i];
          dedx += (*pv1)[j];
        }
        pv->PutValue(j, dedx);
      }
      if(splineFlag) { pv->FillSecondDerivatives(); }
      G4PhysicsTableHelper::SetPhysicsVector(dedxTable, i, pv);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::BuildRangeTable(const G4PhysicsTable* dedxTable,
                                         G4PhysicsTable* rangeTable,
                                         G4bool isIonisation)
// Build range table from the energy loss table
{
  size_t nCouples = dedxTable->size();
  if(0 >= nCouples) { return; }

  size_t n = 100;
  G4double del = 1.0/(G4double)n;

  for (size_t i=0; i<nCouples; ++i) {
    if(isIonisation) {
      if( !(*theFlag)[i] ) { continue; }
    }
    G4PhysicsLogVector* pv = static_cast<G4PhysicsLogVector*>((*dedxTable)[i]);
    size_t npoints = pv->GetVectorLength();
    size_t bin0    = 0;
    G4double elow  = pv->Energy(0);
    G4double ehigh = pv->Energy(npoints-1);
    G4double dedx1 = (*pv)[0];

    //G4cout << "i= " << i << "npoints= " << npoints << " dedx1= " 
    //<< dedx1 << G4endl;

    // protection for specific cases dedx=0
    if(dedx1 == 0.0) {
      for (size_t k=1; k<npoints; ++k) {
        bin0++;
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
    if(npoints < 2) { npoints = 2; }

    delete (*rangeTable)[i];
    G4PhysicsLogVector* v;
    if(0 == bin0) { v = new G4PhysicsLogVector(*pv); }
    else { v = new G4PhysicsLogVector(elow, ehigh, npoints-1); }

    // dedx is exact zero cannot build range table
    if(2 == npoints) {
      v->PutValue(0,1000.);
      v->PutValue(1,2000.);
      G4PhysicsTableHelper::SetPhysicsVector(rangeTable, i, v);
      return;
    }
    v->SetSpline(splineFlag);

    // assumed dedx proportional to beta
    G4double energy1 = v->Energy(0);
    G4double range   = 2.*energy1/dedx1;
    //G4cout << "range0= " << range << G4endl;
    v->PutValue(0,range);

    for (size_t j=1; j<npoints; ++j) {

      G4double energy2 = v->Energy(j);
      G4double de      = (energy2 - energy1) * del;
      G4double energy  = energy2 + de*0.5;
      G4double sum = 0.0;
      //G4cout << "j= " << j << " e1= " << energy1 << " e2= " << energy2 
      //       << " n= " << n << G4endl;
      for (size_t k=0; k<n; ++k) {
        energy -= de;
        dedx1 = pv->Value(energy);
        if(dedx1 > 0.0) { sum += de/dedx1; }
      }
      range += sum;
      v->PutValue(j,range);
      energy1 = energy2;
    }
    if(splineFlag) { v->FillSecondDerivatives(); }
    G4PhysicsTableHelper::SetPhysicsVector(rangeTable, i, v);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LossTableBuilder::BuildInverseRangeTable(const G4PhysicsTable* rangeTable,
                                           G4PhysicsTable* invRangeTable,
                                           G4bool isIonisation)
// Build inverse range table from the energy loss table
{
  size_t nCouples = rangeTable->size();
  if(0 >= nCouples) { return; }

  for (size_t i=0; i<nCouples; ++i) {

    if(isIonisation) {
      if( !(*theFlag)[i] ) { continue; }
    }
    G4PhysicsVector* pv = (*rangeTable)[i];
    size_t npoints = pv->GetVectorLength();
    G4double rlow  = (*pv)[0];
    G4double rhigh = (*pv)[npoints-1];
      
    delete (*invRangeTable)[i];
    G4LPhysicsFreeVector* v = new G4LPhysicsFreeVector(npoints,rlow,rhigh);
    v->SetSpline(splineFlag);

    for (size_t j=0; j<npoints; ++j) {
      G4double e  = pv->Energy(j);
      G4double r  = (*pv)[j];
      v->PutValues(j,r,e);
    }
    if(splineFlag) { v->FillSecondDerivatives(); }

    G4PhysicsTableHelper::SetPhysicsVector(invRangeTable, i, v);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LossTableBuilder::InitialiseBaseMaterials(G4PhysicsTable* table)
{
  size_t nCouples = table->size();
  size_t nFlags = theFlag->size();

  if(nCouples == nFlags && isInitialized) { return; }

  isInitialized = true;

  //G4cout << "%%%%%% G4LossTableBuilder::InitialiseBaseMaterials Ncouples= " 
  //         << nCouples << "  FlagSize= " << nFlags << G4endl; 

  // variable density check
  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();

  /*
  for(size_t i=0; i<nFlags; ++i) {
    G4cout << "CoupleIdx= " << i << "  Flag= " <<  (*theFlag)[i] 
           << " tableFlag= " << table->GetFlag(i) << "  "
           << theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial()->GetName()
           << G4endl;
  }
  */

  // expand vectors
  if(nFlags < nCouples) {
    for(size_t i=nFlags; i<nCouples; ++i) { 
      theDensityFactor->push_back(1.0); 
    }
    for(size_t i=nFlags; i<nCouples; ++i) { theDensityIdx->push_back(-1); }
    for(size_t i=nFlags; i<nCouples; ++i) { theFlag->push_back(true); }
  }
  for(size_t i=0; i<nCouples; ++i) {

    // base material is needed only for a couple which is not
    // initialised and for which tables will be computed
    (*theFlag)[i] = table->GetFlag(i);
    if ((*theDensityIdx)[i] < 0) {
      (*theDensityIdx)[i] = i;
      const G4MaterialCutsCouple* couple =
        theCoupleTable->GetMaterialCutsCouple(i);
      const G4ProductionCuts* pcuts = couple->GetProductionCuts();
      const G4Material* mat  = couple->GetMaterial();
      const G4Material* bmat = mat->GetBaseMaterial();

      // base material exists - find it and check if it can be reused
      if(bmat) {
        for(size_t j=0; j<nCouples; ++j) {

          if(j == i) { continue; }
          const G4MaterialCutsCouple* bcouple =
              theCoupleTable->GetMaterialCutsCouple(j);

          if(bcouple->GetMaterial() == bmat && 
             bcouple->GetProductionCuts() == pcuts) {

            // based couple exist in the same region
            (*theDensityIdx)[i] = j;
            (*theDensityFactor)[i] = mat->GetDensity()/bmat->GetDensity();
            (*theFlag)[i] = false;

            // ensure that there will no double initialisation
            (*theDensityIdx)[j] = j;
            (*theDensityFactor)[j] = 1.0;
            (*theFlag)[j] = true;
            break;
          }
        }
      }
    }
  }
  /*
  for(size_t i=0; i<nCouples; ++i) {
    G4cout << "CoupleIdx= " << i << "  Flag= " <<  (*theFlag)[i] 
           << "  TableFlag= " << table->GetFlag(i) << "  "
           << theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial()->GetName()
           << G4endl;
  }
  G4cout << "%%%%%% G4LossTableBuilder::InitialiseBaseMaterials end" 
         << G4endl; 
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::InitialiseCouples()
{
  isInitialized = true;

  // variable density initialisation for the cas without tables 
  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  size_t nCouples = theCoupleTable->GetTableSize();
  //G4cout << "%%%%%% G4LossTableBuilder::InitialiseCouples() nCouples= " 
  //         << nCouples << G4endl; 

  theDensityFactor->resize(nCouples, 1.0);
  theDensityIdx->resize(nCouples, -1);
  theFlag->resize(nCouples, true);

  for(size_t i=0; i<nCouples; ++i) {

    // base material is needed only for a couple which is not
    // initialised and for which tables will be computed
    if ((*theDensityIdx)[i] < 0) {
      (*theDensityIdx)[i] = i;
      const G4MaterialCutsCouple* couple =
        theCoupleTable->GetMaterialCutsCouple(i);
      const G4ProductionCuts* pcuts = couple->GetProductionCuts();
      const G4Material* mat  = couple->GetMaterial();
      const G4Material* bmat = mat->GetBaseMaterial();

      // base material exists - find it and check if it can be reused
      if(bmat) {
        for(size_t j=0; j<nCouples; ++j) {

          if(j == i) { continue; }
          const G4MaterialCutsCouple* bcouple =
            theCoupleTable->GetMaterialCutsCouple(j);

          if(bcouple->GetMaterial() == bmat && 
             bcouple->GetProductionCuts() == pcuts) {

            // based couple exist in the same region
            (*theDensityIdx)[i] = j;
            (*theDensityFactor)[i] = mat->GetDensity()/bmat->GetDensity();
            (*theFlag)[i] = false;

            // ensure that there will no double initialisation
            (*theDensityIdx)[j] = j;
            (*theDensityFactor)[j] = 1.0;
            (*theFlag)[j] = true;
            break;
          }
        }
      }
    }
  }
  /*  
  for(size_t i=0; i<nCouples; ++i) {
    G4cout << "CoupleIdx= " << i << "  Flag= " <<  (*theFlag)[i] << "  "
           << theCoupleTable->GetMaterialCutsCouple(i)->GetMaterial()->GetName()
           << "  DensityFactor= " << (*theDensityFactor)[i]
           << G4endl;
  }
  G4cout << "%%%%%% G4LossTableBuilder::InitialiseCouples() end" << G4endl; 
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
  if(!table) { return table; }
  if(emin >= emax) { 
    table->clearAndDestroy();
    delete table;
    table = 0;
    return table; 
  }
  InitialiseBaseMaterials(table);

  G4int nbins = theParameters->NumberOfBinsPerDecade();

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  G4PhysicsLogVector* aVector = 0;

  for(size_t i=0; i<numOfCouples; ++i) {

    //G4cout<< "i= " << i << " Flag=  " << GetFlag(i) << G4endl;

    if (GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
        theCoupleTable->GetMaterialCutsCouple(i);
      delete (*table)[i];

      // if start from zero then change the scale

      const G4Material* mat = couple->GetMaterial();

      G4double tmin = std::max(emin,model->MinPrimaryEnergy(mat,part));
      if(0.0 >= tmin) { tmin = eV; }
      G4int n = nbins;

      if(tmin >= emax) {
        aVector = 0;
      } else {
        n *= G4int(std::log10(emax/tmin) + 0.5);
        if(n < 3) { n = 3; }
        aVector = new G4PhysicsLogVector(tmin, emax, n);
      }

      if(aVector) {
        aVector->SetSpline(spline);
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
  return table; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

