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
// File name:     G4LossTableBuilder
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 23.01.2003 V.Ivanchenko Cut per region
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4LossTableBuilder::BuildDEDXTable(
                const G4std::vector<G4PhysicsTable*>& list)    
{   
  size_t n_processes = list.size();

  if      (0 == n_processes) return 0;
  else if (1 == n_processes) return list[0];

  size_t n_vectors = list[0]->length();

  if(!n_vectors) return 0;

  G4PhysicsTable* theTable = new G4PhysicsTable();
  G4bool b;

  for (size_t i=0; i<n_vectors; i++) {

    G4PhysicsVector* pv = (*(list[0]))[i];
    size_t nbins = pv->GetVectorLength();

    G4double elow = pv->GetLowEdgeEnergy(0);
    G4double ehigh = pv->GetLowEdgeEnergy(nbins);
    G4PhysicsLogVector* v = new G4PhysicsLogVector(elow, ehigh, nbins);

    for (size_t j=0; j<nbins; j++) {
      G4double dedx = 0.0;
      G4double energy = pv->GetLowEdgeEnergy(j);

      for (size_t k=0; k<n_processes; k++) {

        dedx += ((*(list[k]))[i])->GetValue(energy, b);

      }
      v->PutValue(j, dedx);
    }
    theTable->insert(v);
  }
  return theTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4LossTableBuilder::BuildRangeTable(
                              const G4PhysicsTable* dedxTable)
// Build range table from the energy loss table
{
  size_t n_vectors = dedxTable->length();
  if(!n_vectors) return 0;

  G4bool b;
  size_t n = 100;
  G4double del = 1.0/(G4double)n;
  G4PhysicsTable* theTable = new G4PhysicsTable();

  for (size_t i=0; i<n_vectors; i++) {

    G4PhysicsVector* pv = (*dedxTable)[i];
    size_t nbins = pv->GetVectorLength();
    G4double elow = pv->GetLowEdgeEnergy(0);
    G4double ehigh = pv->GetLowEdgeEnergy(nbins);
    G4PhysicsLogVector* v = new G4PhysicsLogVector(elow, ehigh, nbins);

    G4double dedx1  = pv->GetValue(elow, b);
    G4double range  = 0.5*elow/dedx1;
    v->PutValue(0,range);
    G4double energy1 = elow;

    for (size_t j=1; j<nbins; j++) {

      G4double energy2 = pv->GetLowEdgeEnergy(j);
      G4double dedx2   = pv->GetValue(energy2, b);
      G4double de = (energy2 - energy1) * del;
      G4double energy = energy1 - de*0.5;
      G4double fac = log(dedx2/dedx1)/log(energy2/energy1);

      for (size_t k=0; k<n; k++) {
        energy += de;
        G4double f = dedx1*exp(fac*log(energy/energy1));
        range  += de/f;
      }
      v->PutValue(j,range);
      energy1 = energy2;
      dedx1   = dedx2;
    }
    theTable->insert(v);
  }

  return theTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4LossTableBuilder::BuildInverseRangeTable(
                      const G4PhysicsTable* dedxTable,
                      const G4PhysicsTable* rangeTable)
// Build inverse range table from the energy loss table
{
  size_t n_vectors = rangeTable->length();
  if(!n_vectors) return 0;

  G4bool b;
  G4PhysicsTable* theTable = new G4PhysicsTable();

  for (size_t i=0; i<n_vectors; i++) {

    G4PhysicsVector* pv = (*rangeTable)[i];
    size_t nbins   = pv->GetVectorLength();
    G4double elow  = pv->GetLowEdgeEnergy(0);
    G4double ehigh = pv->GetLowEdgeEnergy(nbins-1);
    G4double rlow  = pv->GetValue(elow, b);
    G4double rhigh = pv->GetValue(ehigh, b);

    rhigh *= exp(log(rhigh/rlow)/((G4double)(nbins-1)));

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

      G4double e = log(energy1) + log(energy2/energy1)*log(range/range1)/log(range2/range1);

      v->PutValue(j,exp(e));
    }
    theTable->insert(v);
  }

  return theTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
