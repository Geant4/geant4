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
// $Id: G4LossTableBuilder.cc,v 1.27 2008/07/22 15:55:15 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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

  G4bool b;

  G4PhysicsVector* pv = (*(list[0]))[0];
  size_t nbins = pv->GetVectorLength();
  G4double elow = pv->GetLowEdgeEnergy(0);
  G4double ehigh = pv->GetLowEdgeEnergy(nbins);
  for (size_t i=0; i<n_vectors; i++) {

    pv = new G4PhysicsLogVector(elow, ehigh, nbins);
    pv->SetSpline(splineFlag);
    for (size_t j=0; j<nbins; j++) {
      G4double dedx = 0.0;
      G4double energy = pv->GetLowEdgeEnergy(j);

      for (size_t k=0; k<n_processes; k++) {
	dedx += ((*(list[k]))[i])->GetValue(energy, b);
      }
      pv->PutValue(j, dedx);
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
  size_t n_vectors = dedxTable->length();
  if(!n_vectors) return;

  G4bool b;
  size_t n = 100;
  G4double del = 1.0/(G4double)n;

  for (size_t i=0; i<n_vectors; i++) {

    if (rangeTable->GetFlag(i) || !isIonisation) {
      G4PhysicsVector* pv = (*dedxTable)[i];
      size_t nbins = pv->GetVectorLength();
      size_t bin0  = 0;
      G4double elow = pv->GetLowEdgeEnergy(0);
      G4double ehigh = pv->GetLowEdgeEnergy(nbins);
      G4double dedx1  = pv->GetValue(elow, b);

      // protection for specific cases dedx=0
      if(dedx1 == 0.0) {
        for (size_t k=1; k<nbins; k++) {
          bin0++;
          elow  = pv->GetLowEdgeEnergy(k);
          dedx1 = pv->GetValue(elow, b);
          if(dedx1 > 0.0) break;
        }
        nbins -= bin0;
      }

      // initialisation of a new vector
      G4PhysicsLogVector* v = new G4PhysicsLogVector(elow, ehigh, nbins);
      v->SetSpline(splineFlag);

      // assumed dedx proportional to beta
      G4double range  = 2.*elow/dedx1;
      v->PutValue(0,range);
      G4double energy1 = elow;

      for (size_t j=1; j<nbins; j++) {

        G4double energy2 = pv->GetLowEdgeEnergy(j+bin0);
        G4double de      = (energy2 - energy1) * del;
        G4double energy  = energy2 + de*0.5;
        for (size_t k=0; k<n; k++) {
          energy -= de;
          dedx1 = pv->GetValue(energy, b);
          if(dedx1 > 0.0) range += de/dedx1;
	}

	//	G4cout << "Range i= " <<i << " j= " << j << G4endl;
        v->PutValue(j,range);
        energy1 = energy2;
      }
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
  G4bool b;

  for (size_t i=0; i<n_vectors; i++) {

    if (invRangeTable->GetFlag(i) || !isIonisation) {
      G4PhysicsVector* pv = (*rangeTable)[i];
      size_t nbins   = pv->GetVectorLength();
      G4double elow  = pv->GetLowEdgeEnergy(0);
      G4double ehigh = pv->GetLowEdgeEnergy(nbins-1);
      G4double rlow  = pv->GetValue(elow, b);
      G4double rhigh = pv->GetValue(ehigh, b);
      
      G4LPhysicsFreeVector* v = new G4LPhysicsFreeVector(nbins,rlow,rhigh);
      v->SetSpline(splineFlag);

      for (size_t j=0; j<nbins; j++) {
	G4double e  = pv->GetLowEdgeEnergy(j);
	G4double r  = pv->GetValue(e, b);
        v->PutValues(j,r,e);
      }
      v->PutValues(nbins,rhigh+rlow,ehigh);

      G4PhysicsTableHelper::SetPhysicsVector(invRangeTable, i, v);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

