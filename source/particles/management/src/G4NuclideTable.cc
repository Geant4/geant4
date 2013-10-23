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
// MODULE:              G4NuclideTable.cc
//
// Date:                10/10/13
// Author:              T.Koi
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HISTORY
// Based on G4IsomerTable
////////////////////////////////////////////////////////////////////////////////
//
#include "G4NuclideTable.hh"

#include "G4ios.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>
#include <fstream>
#include <sstream>

//const G4double G4NuclideTable::levelTolerance = 1.0*keV;
const G4double G4NuclideTable::levelTolerance = 1.0e-3*eV;
//  torelance for excitation energy
  
 
///////////////////////////////////////////////////////////////////////////////
G4NuclideTable::G4NuclideTable()
  :G4VIsotopeTable("Isomer"),
   threshold_of_half_life(1.0*ns),
   fUserDefinedList(NULL), 
   fIsotopeList(0) 
{
  SetVerboseLevel(G4ParticleTable::GetParticleTable()->GetVerboseLevel());
  //FillIsotopeList();
  return;
}

///////////////////////////////////////////////////////////////////////////////
void G4NuclideTable::FillIsotopeList()
{
  if(fIsotopeList !=0) return;
  
  fIsotopeList = new G4IsotopeList();
  fIsotopeList->reserve(nEntries);

/*
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
*/
}

///////////////////////////////////////////////////////////////////////////////
G4NuclideTable::~G4NuclideTable()
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
//G4NuclideTable::G4NuclideTable(const  G4NuclideTable & right)   
//  :G4VIsotopeTable(right),
//   fIsotopeList(0)
//{
//  FillIsotopeList();
//}

///////////////////////////////////////////////////////////////////////////////
//G4NuclideTable & G4NuclideTable::operator= (const  G4NuclideTable &)
//{
//  return *this;
//}

///////////////////////////////////////////////////////////////////////////////
//
G4IsotopeProperty* G4NuclideTable::GetIsotope(G4int Z, G4int A, G4double E)
{

   G4IsotopeProperty* fProperty = 0;
   G4int ionCode = 1000*Z + A;
   if ( map_pre_load_list.find( ionCode ) !=  map_pre_load_list.end() ) {

      std::multimap< G4double , G4IsotopeProperty* >::iterator lower_bound_itr = 
      map_pre_load_list.find( ionCode ) -> second.lower_bound ( E );

      //std::multimap< G4double , G4IsotopeProperty* >::iterator upper_bound_itr = 
      //map_pre_load_list.find( ionCode ) -> second.upper_bound ( E );

      if ( lower_bound_itr->first - levelTolerance/2 < E && E < lower_bound_itr->first + levelTolerance/2 ) {
         return lower_bound_itr->second;
      }
   }

   char* path = getenv("G4ENSDFSTATEDATA");

   if ( !path ) {
      return fProperty; // not found;
   }

   return fProperty; // not found;

/* 
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
*/
  
}

///////////////////////////////////////////////////////////////////////
G4IsotopeProperty* 
 G4NuclideTable::GetIsotopeByIsoLvl(G4int Z, G4int A, G4int lvl)
{
  if(lvl==0) return GetIsotope(Z,A,0.0);
  G4IsotopeProperty* fProperty = 0;
/*
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
*/
  return fProperty;
  
}

///////////////////////////////////////////////////////////////////////////////
void G4NuclideTable::GenerateNuclide()
{

   if( fIsotopeList !=0 ) return;
  
   fIsotopeList = new G4IsotopeList();

   for (size_t i=0; i<nEntries_ground_state; i++) {

      G4int    ionZ     = groundStateTable[i][idxZ];
      G4int    ionA     = groundStateTable[i][idxA];
      G4int    lvl      = 0; // ground state
      G4double ionE     = groundStateTable[i][idxEnergy]*keV;
      G4double ionLife  = groundStateTable[i][idxLife]*ns;
      G4int    ionJ     = (G4int)(groundStateTable[i][idxSpin]);
      G4double ionMu    = groundStateTable[i][idxMu]*(joule/tesla);

      if ( ionLife == -1 || ionLife > threshold_of_half_life ) {

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

         G4int ionCode = 1000*ionZ + ionA;
         if ( map_pre_load_list.find ( ionCode ) == map_pre_load_list.end() ) {
            std::multimap<G4double, G4IsotopeProperty*> aMultiMap;
            map_pre_load_list.insert( std::pair< G4int , std::multimap< G4double , G4IsotopeProperty* > > ( ionCode , aMultiMap ) );
         }
         map_pre_load_list.find ( ionCode ) -> second.insert( std::pair< G4double, G4IsotopeProperty* >( ionE , fProperty ) );

      }
   }

   if ( threshold_of_half_life >= 1.0*ns ) {

      G4int ionCode=0;
      G4int iLevel=0;
      
      for (size_t i=0; i<nEntries_excite_state; i++) {

         G4int    ionZ     = exciteStateTable[i][idxZ];
         G4int    ionA     = exciteStateTable[i][idxA];
         if ( ionCode != 1000*ionZ + ionA ) {
            iLevel = 0;
            ionCode = 1000*ionZ + ionA;
         } 

         G4double ionE     = exciteStateTable[i][idxEnergy]*keV;
         G4double ionLife  = exciteStateTable[i][idxLife]*ns;
         G4int    ionJ     = (G4int)(exciteStateTable[i][idxSpin]);
         G4double ionMu    = exciteStateTable[i][idxMu]*(joule/tesla);

         if ( ionLife == -1 || ionLife > threshold_of_half_life ) {

            iLevel++;
            if ( iLevel > 9 ) iLevel=9;
            //G4cout << ionZ << " " << ionA << " " << iLevel << " " << ionE/keV << " [keV]" << G4endl;

            G4IsotopeProperty* fProperty = new G4IsotopeProperty(); 

            // Set Isotope Property
            fProperty->SetAtomicNumber(ionZ);
            fProperty->SetAtomicMass(ionA);
            fProperty->SetIsomerLevel(iLevel);
            fProperty->SetEnergy(ionE);
            fProperty->SetiSpin(ionJ);
            fProperty->SetLifeTime(ionLife);
            fProperty->SetDecayTable(0);
            fProperty->SetMagneticMoment(ionMu);
       
            fIsotopeList->push_back(fProperty);

            if ( map_pre_load_list.find ( ionCode ) == map_pre_load_list.end() ) {
               std::multimap<G4double, G4IsotopeProperty*> aMultiMap;
               map_pre_load_list.insert( std::pair< G4int , std::multimap< G4double , G4IsotopeProperty* > > ( ionCode , aMultiMap ) );
            }
            map_pre_load_list.find ( ionCode ) -> second.insert( std::pair< G4double, G4IsotopeProperty* >( ionE , fProperty ) );

         }
      }
   } else {

      char* path = getenv("G4ENSDFSTATEDATA");

      if ( !path ) {
         G4Exception("G4NuclideTable", "XXX",
                  FatalException, "G4ENSDFSTATEDATA environment variable must be set");
      }
   
      std::fstream ifs;
      G4String filename(path);
      filename += "/ENSDFSTATE.dat";

      ifs.open( filename.c_str() );
     
      if ( !ifs.good() ) {
         G4Exception("G4NuclideTable", "XXX",
                  FatalException, "Not find ENSDFSTATE.dat");
      }
     

      G4int ionCode=0;
      G4int iLevel=0;

      G4int ionZ;
      G4int ionA;
      G4double ionE;
      G4double ionLife;
      G4int ionJ;
      G4double ionMu;
      
      ifs >> ionZ >> ionA >> ionE >> ionLife >> ionJ >> ionMu;

      while ( ifs.good() ) {

         if ( ionCode != 1000*ionZ + ionA ) {
            iLevel = 0;
            ionCode = 1000*ionZ + ionA;
         } 

         ionE *= keV;
         ionLife *= ns;
         ionMu *= (joule/tesla);

         if ( ionLife == -1 || ionLife > threshold_of_half_life ) {

            iLevel++;
            if ( iLevel > 9 ) iLevel=9;
            G4cout << ionZ << " " << ionA << " " << iLevel << " " << ionE/keV << " [keV]" << G4endl;

            G4IsotopeProperty* fProperty = new G4IsotopeProperty(); 

            // Set Isotope Property
            fProperty->SetAtomicNumber(ionZ);
            fProperty->SetAtomicMass(ionA);
            fProperty->SetIsomerLevel(iLevel);
            fProperty->SetEnergy(ionE);
            fProperty->SetiSpin(ionJ);
            fProperty->SetLifeTime(ionLife);
            fProperty->SetDecayTable(0);
            fProperty->SetMagneticMoment(ionMu);
       
            fIsotopeList->push_back(fProperty);

            if ( map_pre_load_list.find ( ionCode ) == map_pre_load_list.end() ) {
               std::multimap<G4double, G4IsotopeProperty*> aMultiMap;
               map_pre_load_list.insert( std::pair< G4int , std::multimap< G4double , G4IsotopeProperty* > > ( ionCode , aMultiMap ) );
            }
            map_pre_load_list.find ( ionCode ) -> second.insert( std::pair< G4double, G4IsotopeProperty* >( ionE , fProperty ) );

         }

         ifs >> ionZ >> ionA >> ionE >> ionLife >> ionJ >> ionMu;
      }

   }

   if ( fUserDefinedList != NULL ) {
      for ( G4IsotopeList::iterator it = fUserDefinedList->begin() ; it != fUserDefinedList->end() ; it++ ) {
         fIsotopeList->push_back( *it );
      }
   }

}


void G4NuclideTable::AddState( G4int ionZ, G4int ionA, G4double ionE, G4double ionLife, G4int ionJ=0, G4double ionMu=0.0)
{
   if ( fUserDefinedList == NULL ) fUserDefinedList = new G4IsotopeList();

            G4IsotopeProperty* fProperty = new G4IsotopeProperty(); 

            // Set Isotope Property
            fProperty->SetAtomicNumber(ionZ);
            fProperty->SetAtomicMass(ionA);
            fProperty->SetIsomerLevel(9);
            fProperty->SetEnergy(ionE);
            fProperty->SetiSpin(ionJ);
            fProperty->SetLifeTime(ionLife);
            fProperty->SetDecayTable(0);
            fProperty->SetMagneticMoment(ionMu);
       
            fUserDefinedList->push_back(fProperty);

}

