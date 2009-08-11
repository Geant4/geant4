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
// $Id: G4LossTableBuilder.cc,v 1.32 2009-08-11 17:24:53 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4LPhysicsFreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LossTableBuilder::G4LossTableBuilder() 
{
  splineFlag = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LossTableBuilder::~G4LossTableBuilder() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LossTableBuilder::BuildDEDXTable(G4PhysicsTable* dedxTable,
				   const std::vector<G4PhysicsTable*>& list)
{
  size_t n_processes = list.size();
  if(1 >= n_processes) return;

  size_t n_vectors = (list[0])->length();
  if(0 >= n_vectors) return;

  G4PhysicsLogVector* pv0 = static_cast<G4PhysicsLogVector*>((*(list[0]))[0]);
  size_t npoints = pv0->GetVectorLength();
  for (size_t i=0; i<n_vectors; i++) {

    G4PhysicsLogVector* pv = new G4PhysicsLogVector(*pv0);
    //    pv = new G4PhysicsLogVector(elow, ehigh, npoints-1);
    pv->SetSpline(splineFlag);
    for (size_t j=0; j<npoints; j++) {
      G4double dedx = 0.0;
      for (size_t k=0; k<n_processes; k++) {
        G4PhysicsVector* pv1   = (*(list[k]))[i];
	dedx += (*pv1)[j];
      }
      pv->PutValue(j, dedx);
    }
    if(splineFlag) pv->FillSecondDerivatives();
    G4PhysicsTableHelper::SetPhysicsVector(dedxTable, i, pv);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::BuildRangeTable(const G4PhysicsTable* dedxTable,
					 G4PhysicsTable* rangeTable,
					 G4bool isIonisation)
// Build range table from the energy loss table
{
  size_t n_vectors = dedxTable->length();
  if(!n_vectors) return;

  size_t n = 100;
  G4double del = 1.0/(G4double)n;

  for (size_t i=0; i<n_vectors; i++) {

    if (rangeTable->GetFlag(i) || !isIonisation) {
      G4PhysicsLogVector* pv = 
	static_cast<G4PhysicsLogVector*>((*dedxTable)[i]);
      size_t npoints = pv->GetVectorLength();
      size_t bin0    = 0;
      G4double elow  = pv->Energy(0);
      G4double ehigh = pv->Energy(npoints-1);
      G4double dedx1 = pv->Value(elow);

      //G4cout << "nbins= " << nbins << " dedx1= " << dedx1 << G4endl;

      // protection for specific cases dedx=0
      if(dedx1 == 0.0) {
        for (size_t k=1; k<npoints; k++) {
          bin0++;
          elow  = pv->Energy(k);
          dedx1 = (*pv)[k];
          if(dedx1 > 0.0) break;
        }
        npoints -= bin0;
      }
      // G4cout<<"New Range vector" << G4endl;
      // G4cout<<"nbins= "<<npoints-1<<" elow= "<<elow<<" ehigh= "<<ehigh<<G4endl;
      // initialisation of a new vector
      if(npoints < 2) npoints = 2;
      G4PhysicsLogVector* v;
      if(0 == bin0) { v = new G4PhysicsLogVector(*pv); }
      else          { v = new G4PhysicsLogVector(elow, ehigh, npoints-1); }
      // dedx is exect zero
      if(2 == npoints) {
	v->PutValue(0,1000.);
	v->PutValue(1,2000.);
	G4PhysicsTableHelper::SetPhysicsVector(rangeTable, i, v);
	return;
      }
      v->SetSpline(splineFlag);

      // assumed dedx proportional to beta
      G4double range  = 2.*elow/dedx1;
      v->PutValue(0,range);
      G4double energy1 = elow;

      for (size_t j=1; j<npoints; j++) {

        G4double energy2 = pv->Energy(j+bin0);
        G4double de      = (energy2 - energy1) * del;
        G4double energy  = energy2 + de*0.5;
        G4double sum = 0.0;
        for (size_t k=0; k<n; k++) {
          energy -= de;
          dedx1 = pv->Value(energy);
          if(dedx1 > 0.0) sum += de/dedx1;
	}
        range += sum;
	//	G4cout << "Range i= " <<i << " j= " << j << G4endl;
        v->PutValue(j,range);
        energy1 = energy2;
      }
      if(splineFlag) v->FillSecondDerivatives();
      G4PhysicsTableHelper::SetPhysicsVector(rangeTable, i, v);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::BuildInverseRangeTable(const G4PhysicsTable* rangeTable,
						G4PhysicsTable* invRangeTable,
						G4bool isIonisation)
// Build inverse range table from the energy loss table
{
  size_t n_vectors = rangeTable->length();
  if(!n_vectors) return;

  for (size_t i=0; i<n_vectors; i++) {

    if (invRangeTable->GetFlag(i) || !isIonisation) {
      G4PhysicsVector* pv = (*rangeTable)[i];
      size_t npoints = pv->GetVectorLength();
      G4double rlow  = (*pv)[0];
      G4double rhigh = (*pv)[npoints-1];
      
      G4LPhysicsFreeVector* v = new G4LPhysicsFreeVector(npoints,rlow,rhigh);
      v->SetSpline(splineFlag);

      for (size_t j=0; j<npoints; j++) {
	G4double e  = pv->Energy(j);
	G4double r  = (*pv)[j];
        v->PutValues(j,r,e);
      }
      if(splineFlag) v->FillSecondDerivatives();

      G4PhysicsTableHelper::SetPhysicsVector(invRangeTable, i, v);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

