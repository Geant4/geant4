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
// $Id: G4LossTableBuilder.cc,v 1.17 2005/04/12 18:13:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 07-12-04 Fix of BuildDEDX table (V.Ivantchenko)
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::BuildDEDXTable(G4PhysicsTable* dedxTable,
                                  const std::vector<G4PhysicsTable*>& list)
{
  size_t n_processes = list.size();
  if(1 >= n_processes) return;

  size_t n_vectors = dedxTable->length();
  if(0 >= n_vectors) return;

  G4bool b;

  for (size_t i=0; i<n_vectors; i++) {

    G4PhysicsVector* pv = (*dedxTable)[i];
    size_t nbins = pv->GetVectorLength();

    for (size_t j=0; j<nbins; j++) {
      G4double dedx = 0.0;
      G4double energy = pv->GetLowEdgeEnergy(j);

      for (size_t k=0; k<n_processes; k++) {
	dedx += ((*(list[k]))[i])->GetValue(energy, b);
      }
      pv->PutValue(j, dedx);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::BuildRangeTable(const G4PhysicsTable* dedxTable,
                                               G4PhysicsTable* rangeTable)
// Build range table from the energy loss table
{
  size_t n_vectors = dedxTable->length();
  if(!n_vectors) return;

  G4bool b;
  size_t n = 100;
  G4double del = 1.0/(G4double)n;

  for (size_t i=0; i<n_vectors; i++) {

    if (rangeTable->GetFlag(i)) {
      G4PhysicsVector* pv = (*dedxTable)[i];
      size_t nbins = pv->GetVectorLength();
      size_t bin0  = 0;
      G4double elow = pv->GetLowEdgeEnergy(0);
      G4double ehigh = pv->GetLowEdgeEnergy(nbins);
      G4double dedx1  = pv->GetValue(elow, b);

      if(dedx1 == 0.0) {
        for (size_t k=1; k<nbins; k++) {
          bin0++;
          elow  = pv->GetLowEdgeEnergy(k);
          dedx1 = pv->GetValue(elow, b);
          if(dedx1 > 0.0) break;
        }
        nbins -= bin0;
      }

      G4PhysicsLogVector* v = new G4PhysicsLogVector(elow, ehigh, nbins);

      G4double range  = 2.*elow/dedx1;
      //G4double range  = elow/dedx1;
      v->PutValue(0,range);
      G4double energy1 = elow;

      for (size_t j=1; j<nbins; j++) {

        G4double energy2 = pv->GetLowEdgeEnergy(j+bin0);
        G4double dedx2   = pv->GetValue(energy2, b);
        G4double de      = (energy2 - energy1) * del;
        G4double energy  = energy1 - de*0.5;

        G4bool   yes     = true;
        if(dedx1 < DBL_MIN || dedx2 < DBL_MIN) yes = false;   

        G4double fac, f;

        if(yes) fac = std::log(dedx2/dedx1)/std::log(energy2/energy1);
        else    fac = (dedx2 - dedx1)/(energy2 - energy1);

        for (size_t k=0; k<n; k++) {
          energy += de;
          if(yes) f = dedx1*std::exp(fac*std::log(energy/energy1));
          else    f = dedx1 + fac*(energy - energy1);
          if(f > DBL_MIN) range  += de/f;
        }
        v->PutValue(j,range);
        energy1 = energy2;
        dedx1   = dedx2;
      }
      G4PhysicsTableHelper::SetPhysicsVector(rangeTable, i, v);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableBuilder::BuildInverseRangeTable(const G4PhysicsTable* rangeTable,
                                                      G4PhysicsTable* invRangeTable)
// Build inverse range table from the energy loss table
{
  size_t n_vectors = rangeTable->length();
  if(!n_vectors) return;
  G4bool b;

  for (size_t i=0; i<n_vectors; i++) {

    if (invRangeTable->GetFlag(i)) {
      G4PhysicsVector* pv = (*rangeTable)[i];
      size_t nbins   = pv->GetVectorLength();
      G4double elow  = pv->GetLowEdgeEnergy(0);
      G4double ehigh = pv->GetLowEdgeEnergy(nbins-1);
      G4double rlow  = pv->GetValue(elow, b);
      G4double rhigh = pv->GetValue(ehigh, b);

      rhigh *= std::exp(std::log(rhigh/rlow)/((G4double)(nbins-1)));

      G4PhysicsLogVector* v = new G4PhysicsLogVector(rlow, rhigh, nbins);

      v->PutValue(0,elow);
      G4double energy1 = elow;
      G4double range1  = rlow;
      G4double energy2 = elow;
      G4double range2  = rlow;
      size_t ilow      = 0;
      size_t ihigh;

      for (size_t j=1; j<nbins; j++) {

        G4double range = v->GetLowEdgeEnergy(j);

        for (ihigh=ilow+1; ihigh<nbins; ihigh++) {
          energy2 = pv->GetLowEdgeEnergy(ihigh);
          range2  = pv->GetValue(energy2, b);
          if(range2 >= range || ihigh == nbins-1) {
            ilow = ihigh - 1;
            energy1 = pv->GetLowEdgeEnergy(ilow);
            range1  = pv->GetValue(energy1, b);
            break;
	  }
        }

        G4double e = std::log(energy1) + 
                     std::log(energy2/energy1)*std::log(range/range1)/std::log(range2/range1);

        v->PutValue(j,std::exp(e));
      }
      G4PhysicsTableHelper::SetPhysicsVector(invRangeTable, i, v);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

