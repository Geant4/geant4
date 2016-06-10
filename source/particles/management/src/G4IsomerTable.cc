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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4IsomerTable.cc
//
// Date:                5/05/13
// Author:              H.Kurashige
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HISTORY
////////////////////////////////////////////////////////////////////////////////
//
#include "G4IsomerTable.hh"

#include "G4ios.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>
#include <fstream>
#include <sstream>

const G4double G4IsomerTable::levelTolerance = 1.0*keV;
//  torelance for excitation energy
  
 
///////////////////////////////////////////////////////////////////////////////
G4IsomerTable::G4IsomerTable()
  :G4VIsotopeTable("Isomer"),
   fIsotopeList(0)
{
  SetVerboseLevel(G4ParticleTable::GetParticleTable()->GetVerboseLevel());
  FillIsotopeList();
  return;
}

///////////////////////////////////////////////////////////////////////////////
void G4IsomerTable::FillIsotopeList()
{
  if(fIsotopeList !=0) return;
  
  fIsotopeList = new G4IsotopeList();
  fIsotopeList->reserve(nEntries);

  for (size_t i=0; i<nEntries; i++) {
    G4int    pid      = (G4int)(isomerTable[i][idxPID]);
    G4int    ionZ     = pid/10000; 
    G4int    ionA     = (pid-ionZ*10000)/10;
    G4int    lvl      = pid%10;
    G4double ionE     = isomerTable[i][idxEnergy]*keV;
    G4double ionLife  = isomerTable[i][idxLife]*ns;
    G4int    ionJ     = (G4int)(isomerTable[i][idxSpin]);
    G4double ionMu    = isomerTable[i][idxMu]*(joule/tesla);

    G4IsotopeProperty* fProperty = new G4IsotopeProperty(); 
    // Set Isotope Property
    fProperty->SetAtomicNumber(ionZ);
    fProperty->SetAtomicMass(ionA);
    fProperty->SetIsomerLevel(lvl);
    fProperty->SetEnergy(ionE);
    fProperty->SetiSpin(ionJ);
    fProperty->SetLifeTime(ionLife);
    fProperty->SetDecayTable(0);
    fProperty->SetMagneticMoment(ionMu);
    
    fIsotopeList->push_back(fProperty);
  }
}

///////////////////////////////////////////////////////////////////////////////
G4IsomerTable::~G4IsomerTable()
{
  if (fIsotopeList!=0) {
    for (size_t i = 0 ; i<fIsotopeList->size(); i++) {
      delete (*fIsotopeList)[i];
    }
    fIsotopeList->clear();
    delete fIsotopeList;
    fIsotopeList = 0;
  }
}

///////////////////////////////////////////////////////////////////////////////
G4IsomerTable::G4IsomerTable(const  G4IsomerTable & right)   
  :G4VIsotopeTable(right),
   fIsotopeList(0)
{
  FillIsotopeList();
}

///////////////////////////////////////////////////////////////////////////////
G4IsomerTable & G4IsomerTable::operator= (const  G4IsomerTable &)
{
  FillIsotopeList();
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
G4IsotopeProperty* G4IsomerTable::GetIsotope(G4int Z, G4int A, G4double E)
{
  G4IsotopeProperty* fProperty = 0;
  if ((Z<MinZ)||(Z>MaxZ)||(A>MaxA)) return fProperty; // not found
  
  G4int low  = 0;
  G4int high = entries() -1;
  G4int ptr = (low+high)/2;
 
  // find Z
  G4int zptr   = (*fIsotopeList)[ptr]->GetAtomicNumber();
  while ( high-low > 1){
    if (zptr == Z){
      while ((zptr == Z)&&(ptr>0)) {
	ptr -= 1;
	zptr = (*fIsotopeList)[ptr]->GetAtomicNumber();
      }
      if (ptr!=0) ptr += 1;
      break;
    } else if (zptr >Z) {
      high = ptr;
      ptr = (low+high)/2;
    } else {
      low = ptr;
      ptr = (low+high)/2;
    }
    zptr   = (*fIsotopeList)[ptr]->GetAtomicNumber();
  }
  if ( Z == (*fIsotopeList)[low]->GetAtomicNumber()) ptr =low;
  else if ( Z != zptr ) return fProperty; // not found

  // find A
  G4int aptr   = (*fIsotopeList)[ptr]->GetAtomicMass();
  while ( (aptr<A) && (ptr+1<(G4int)(entries())) ) {
    ptr +=1;
    aptr   = (*fIsotopeList)[ptr]->GetAtomicMass();
    if (Z != (*fIsotopeList)[ptr]->GetAtomicNumber()) break;
  }
  if ( aptr != A)  return fProperty;  // not found

  // find E 
  G4double ptrE = (*fIsotopeList)[ptr]->GetEnergy();
  while ( (ptrE < E-levelTolerance ) && (ptr+1<(G4int)(entries())) ){
    ptr +=1;
    ptrE = (*fIsotopeList)[ptr]->GetEnergy();
    if (A != (*fIsotopeList)[ptr]->GetAtomicMass())  return fProperty; // not found
  }
  if (ptrE > E+levelTolerance) return fProperty; // not found
  
  // FOUND !!
  fProperty =  (*fIsotopeList)[ptr];

  return fProperty;
  
}

///////////////////////////////////////////////////////////////////////
G4IsotopeProperty* 
 G4IsomerTable::GetIsotopeByIsoLvl(G4int Z, G4int A, G4int lvl)
{
  G4IsotopeProperty* fProperty = 0;
  if ((Z<MinZ)||(Z>MaxZ)||(A>MaxA)) return fProperty; // not found

  G4int low  = 0;
  G4int high = entries() -1;
  G4int ptr = (low+high)/2;


  // find PID
  const G4int PID = 10000*Z +10*A + lvl; 
  G4int id   =  (fIsotopeList->at(ptr))->GetAtomicNumber() *10000
              + (fIsotopeList->at(ptr))->GetAtomicMass() *10
              + (fIsotopeList->at(ptr))->GetIsomerLevel();
  while ( high-low > 1){
    if (id == PID){
      break;
    } else if (id >PID) {
      high = ptr;
      ptr = (low+high)/2;
    } else {
      low = ptr;
      ptr = (low+high)/2;
    }
    id   =  (fIsotopeList->at(ptr))->GetAtomicNumber() *10000
          + (fIsotopeList->at(ptr))->GetAtomicMass() *10
          + (fIsotopeList->at(ptr))->GetIsomerLevel();
  }

  if (id == PID) {
    fProperty =  (*fIsotopeList)[ptr];
  } else {
    id   =  (fIsotopeList->at(low))->GetAtomicNumber() *10000
          + (fIsotopeList->at(low))->GetAtomicMass() *10
          + (fIsotopeList->at(low))->GetIsomerLevel();
    if (id == PID) {
      fProperty =  (*fIsotopeList)[low];
    } else {
      id   =  (fIsotopeList->at(high))->GetAtomicNumber() *10000
            + (fIsotopeList->at(high))->GetAtomicMass() *10
            + (fIsotopeList->at(high))->GetIsomerLevel();
      if (id == PID) fProperty =  (*fIsotopeList)[high];
    }
  }
  return fProperty;
  
}
