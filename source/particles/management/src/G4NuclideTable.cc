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

const G4double G4NuclideTable::levelTolerance = 1.0*eV;
// const G4double G4NuclideTable::levelTolerance = 1.0e-3*eV;
//  torelance for excitation energy
  
 
///////////////////////////////////////////////////////////////////////////////
G4NuclideTable* G4NuclideTable::GetInstance() {
   static G4NuclideTable instance;
   return &instance;
}

///////////////////////////////////////////////////////////////////////////////
G4NuclideTable::G4NuclideTable()
  :G4VIsotopeTable("Isomer"),
   threshold_of_half_life(1000.0*ns),
   fUserDefinedList(NULL), 
   fIsotopeList(0) 
{
  //SetVerboseLevel(G4ParticleTable::GetParticleTable()->GetVerboseLevel());
  FillHardCodeList();
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

  for ( std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > >::iterator 
     it = map_pre_load_list.begin(); it != map_pre_load_list.end(); it++ ) {
     it->second.clear();
  }
  map_pre_load_list.clear();

  for ( std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > >::iterator 
     it = map_hard_code_list.begin(); it != map_hard_code_list.end(); it++ ) {
     for ( std::multimap< G4double , G4IsotopeProperty* >::iterator 
        itt = it->second.begin(); itt != it->second.end(); itt++ ) {
        delete itt->second;
     }
     it->second.clear();
  }
  map_hard_code_list.clear();

  for ( std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > >::iterator 
     it = map_full_list.begin(); it != map_full_list.end(); it++ ) {
     for ( std::multimap< G4double , G4IsotopeProperty* >::iterator 
        itt = it->second.begin(); itt != it->second.end(); itt++ ) {
        delete itt->second;
     }
     it->second.clear();
  }
  map_full_list.clear();
}

///////////////////////////////////////////////////////////////////////////////
//
G4IsotopeProperty* G4NuclideTable::GetIsotope(G4int Z, G4int A, G4double E)
{

   G4IsotopeProperty* fProperty = 0;
   G4int ionCode = 1000*Z + A;

   //Serching pre-load
   //Note: isomer level is properly set only for pre_load_list.
   if ( map_pre_load_list.find( ionCode ) !=  map_pre_load_list.end() ) {

     std::multimap< G4double , G4IsotopeProperty* >::iterator lower_bound_itr = 
       map_pre_load_list.find( ionCode ) -> second.lower_bound ( E - levelTolerance/2 );
     
     //std::multimap< G4double , G4IsotopeProperty* >::iterator upper_bound_itr = 
     //map_pre_load_list.find( ionCode ) -> second.upper_bound ( E );
     
     G4double levelE = DBL_MAX;
     if ( lower_bound_itr !=  map_pre_load_list.find( ionCode ) -> second.end() ) {
       levelE = lower_bound_itr->first;
       if ( levelE - levelTolerance/2 <= E && E < levelE + levelTolerance/2 ) {
         return lower_bound_itr->second; // found
       } 
     }
   }
     
   //Searching hard-code
   if ( map_hard_code_list.find( ionCode ) !=  map_hard_code_list.end() ) {
      std::multimap< G4double , G4IsotopeProperty* >::iterator lower_bound_itr = 
      map_hard_code_list.find( ionCode ) -> second.lower_bound ( E - levelTolerance/2 );

      //std::multimap< G4double , G4IsotopeProperty* >::iterator upper_bound_itr = 
      //map_pre_load_list.find( ionCode ) -> second.upper_bound ( E );
      
      G4double levelE = DBL_MAX;
      if ( lower_bound_itr !=  map_hard_code_list.find( ionCode ) -> second.end() ) {
         levelE = lower_bound_itr->first;
	 if ( levelE - levelTolerance/2 <= E && E < levelE + levelTolerance/2 ) {
	   return lower_bound_itr->second; // found
	 }
      }
   }

   //Searching big-list 
   char* path = getenv("G4ENSDFSTATEDATA");

   if ( !path ) {
      return fProperty; // not found;
   }

   if ( map_full_list.find( ionCode ) ==  map_full_list.end() ) {

      std::multimap<G4double, G4IsotopeProperty*> aMultiMap;
      map_full_list.insert( std::pair< G4int , std::multimap< G4double , G4IsotopeProperty* > > ( ionCode , aMultiMap ) );

      std::fstream ifs;
      G4String filename(path);
      filename += "/ENSDFSTATE.dat";
      ifs.open( filename.c_str() );

      G4bool reading_target = false; 

      G4int ionZ;
      G4int ionA;
      G4double ionE;
      G4double ionLife;
      G4int ionJ;
      G4double ionMu;
      
      ifs >> ionZ >> ionA >> ionE >> ionLife >> ionJ >> ionMu;

      while ( ifs.good() ) {

         if ( ionZ == Z && ionA == A ) {

            reading_target = true;

            ionE *= keV;
            ionLife *= ns;
            ionMu *= (joule/tesla);

            G4IsotopeProperty* property = new G4IsotopeProperty(); 

            G4int iLevel=9;
            property->SetAtomicNumber(ionZ);
            property->SetAtomicMass(ionA);
            property->SetIsomerLevel(iLevel);
            property->SetEnergy(ionE);
            property->SetiSpin(ionJ);
            property->SetLifeTime(ionLife);
            property->SetMagneticMoment(ionMu);
       
            map_full_list.find ( ionCode ) -> second.insert( std::pair< G4double, G4IsotopeProperty* >( ionE , property ) );

         } else if ( reading_target == true ) {
            ifs.close();
            break;
         }
         
         ifs >> ionZ >> ionA >> ionE >> ionLife >> ionJ >> ionMu;
      }

      ifs.close();
   }


   if ( map_full_list.find( ionCode ) !=  map_full_list.end() ) {

      std::multimap< G4double , G4IsotopeProperty* >::iterator lower_bound_itr = 
      map_full_list.find( ionCode ) -> second.lower_bound ( E - levelTolerance/2 );

      //std::multimap< G4double , G4IsotopeProperty* >::iterator upper_bound_itr = 
      //map_full_list.find( ionCode ) -> second.upper_bound ( E - levelTolerance/2 );
      
      G4double levelE = DBL_MAX;
      if ( lower_bound_itr !=  map_full_list.find( ionCode ) -> second.end() ) {
         levelE = lower_bound_itr->first;
	 if ( levelE - levelTolerance/2 < E && E < levelE + levelTolerance/2 ) {
	   return lower_bound_itr->second; // found
	 }
      }
   }

   return fProperty; // not found;
}

///////////////////////////////////////////////////////////////////////
G4IsotopeProperty* 
 G4NuclideTable::GetIsotopeByIsoLvl(G4int Z, G4int A, G4int lvl)
{
  if(lvl==0) return GetIsotope(Z,A,0.0);
  return (G4IsotopeProperty*)0;
}

///////////////////////////////////////////////////////////////////////////////
void G4NuclideTable::FillHardCodeList()
{
   for (size_t i=0; i<nEntries_ground_state; i++) {

      G4int    ionZ     = (G4int)groundStateTable[i][idxZ];
      G4int    ionA     = (G4int)groundStateTable[i][idxA];
      G4int    lvl      = 0; // ground state
      G4double ionE     = groundStateTable[i][idxEnergy]*keV;
      G4double ionLife  = groundStateTable[i][idxLife]*ns;
      G4int    ionJ     = (G4int)(groundStateTable[i][idxSpin]);
      G4double ionMu    = groundStateTable[i][idxMu]*(joule/tesla);

      G4int ionCode = 1000*ionZ + ionA;

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

      if ( map_hard_code_list.find ( ionCode ) == map_hard_code_list.end() ) {
         std::multimap<G4double, G4IsotopeProperty*> aMultiMap;
         map_hard_code_list.insert( std::pair< G4int , std::multimap< G4double , G4IsotopeProperty* > > ( ionCode , aMultiMap ) );
      }
      map_hard_code_list.find ( ionCode ) -> second.insert( std::pair< G4double, G4IsotopeProperty* >( ionE , fProperty ) );

   }

   for (size_t i=0; i<nEntries_excite_state; i++) {

      G4int    ionZ     = (G4int)exciteStateTable[i][idxZ];
      G4int    ionA     = (G4int)exciteStateTable[i][idxA];
      G4double ionE     = exciteStateTable[i][idxEnergy]*keV;
      G4double ionLife  = exciteStateTable[i][idxLife]*ns;
      G4int    ionJ     = (G4int)(exciteStateTable[i][idxSpin]);
      G4double ionMu    = exciteStateTable[i][idxMu]*(joule/tesla);

      G4int ionCode = 1000*ionZ + ionA;

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

      if ( map_hard_code_list.find ( ionCode ) == map_hard_code_list.end() ) {
         std::multimap<G4double, G4IsotopeProperty*> aMultiMap;
         map_hard_code_list.insert( std::pair< G4int , std::multimap< G4double , G4IsotopeProperty* > > ( ionCode , aMultiMap ) );
      }
      map_hard_code_list.find ( ionCode ) -> second.insert( std::pair< G4double, G4IsotopeProperty* >( ionE , fProperty ) );

   }
}

///////////////////////////////////////////////////////////////////////////////
void G4NuclideTable::GenerateNuclide()
{

   if( fIsotopeList !=0 ) return;
   fIsotopeList = new G4IsotopeList();

   for (size_t i=0; i<nEntries_ground_state; i++) {

      G4int    ionZ     = (G4int)groundStateTable[i][idxZ];
      G4int    ionA     = (G4int)groundStateTable[i][idxA];
      G4int    lvl      = 0; // ground state
      G4double ionE     = groundStateTable[i][idxEnergy]*keV;
      G4double ionLife  = groundStateTable[i][idxLife]*ns;
      G4int    ionJ     = (G4int)(groundStateTable[i][idxSpin]);
      G4double ionMu    = groundStateTable[i][idxMu]*(joule/tesla);

      if ( ionLife < 0.0 || ionLife*std::log(2.0) > threshold_of_half_life ) {

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
    
         //G4cout << ionZ << " " << ionA << " " << lvl << " " << ionE/keV << " [keV]" << G4endl;
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
      G4double previousE=0.0;
      
      for (size_t i=0; i<nEntries_excite_state; i++) {

         G4int    ionZ     = (G4int)exciteStateTable[i][idxZ];
         G4int    ionA     = (G4int)exciteStateTable[i][idxA];
         if ( ionCode != 1000*ionZ + ionA ) {
            previousE=0.0;
            iLevel = 0;
            ionCode = 1000*ionZ + ionA;
         } 

         G4double ionE     = exciteStateTable[i][idxEnergy]*keV;
         G4double ionLife  = exciteStateTable[i][idxLife]*ns;
         G4int    ionJ     = (G4int)(exciteStateTable[i][idxSpin]);
         G4double ionMu    = exciteStateTable[i][idxMu]*(joule/tesla);

         if (( ionLife < 0.0 || ionLife*ionLife*std::log(2.0) > threshold_of_half_life )
           && (ionE > levelTolerance+previousE)) {
            previousE = ionE;
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
         G4Exception("G4NuclideTable", "PART70000",
                  FatalException, "G4ENSDFSTATEDATA environment variable must be set");
	 return;
      }
   
      std::fstream ifs;
      G4String filename(path);
      filename += "/ENSDFSTATE.dat";

      ifs.open( filename.c_str() );
     
      if ( !ifs.good() ) {
         G4Exception("G4NuclideTable", "PART70001",
                  FatalException, "ENSDFSTATE.dat is not found.");
	 return;
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

         //if ( ionLife == -1 || ionLife > threshold_of_half_life ) {
         if ( ionLife*std::log(2.0) > threshold_of_half_life && ionE != 0 ) {

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

