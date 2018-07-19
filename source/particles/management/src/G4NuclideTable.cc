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
#include "G4NuclideTableMessenger.hh"

#include "G4ios.hh"
#include "G4String.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>
#include <fstream>
#include <sstream>

// const G4double G4NuclideTable::levelTolerance = 
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
   minimum_threshold_of_half_life(DBL_MAX),
   fUserDefinedList(NULL), 
   fIsotopeList(NULL),
   flevelTolerance(1.0*eV)
{
  //SetVerboseLevel(G4ParticleTable::GetParticleTable()->GetVerboseLevel());
  //FillHardCodeList();
  fMessenger = new G4NuclideTableMessenger(this);
  fIsotopeList = new G4IsotopeList();
  GenerateNuclide();
}

///////////////////////////////////////////////////////////////////////////////
G4NuclideTable::~G4NuclideTable()
{

  for ( std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > >::iterator 
     it = map_pre_load_list.begin(); it != map_pre_load_list.end(); it++ ) {
     it->second.clear();
  }
  map_pre_load_list.clear();


  for ( std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > >::iterator 
     it = map_full_list.begin(); it != map_full_list.end(); it++ ) {
     it->second.clear();
  }
  map_full_list.clear();

  if (fIsotopeList!=0) {
    for (size_t i = 0 ; i<fIsotopeList->size(); i++) {
       //G4IsotopeProperty* fProperty = (*fIsotopeList)[i]; std::cout << fProperty->GetAtomicNumber() << " " << fProperty->GetAtomicMass() << " " << fProperty->GetEnergy() << std::endl;
      delete (*fIsotopeList)[i];
    }
    fIsotopeList->clear();
    delete fIsotopeList;
    fIsotopeList = 0;
  }

}

///////////////////////////////////////////////////////////////////////////////
//
G4IsotopeProperty* G4NuclideTable::GetIsotope(G4int Z, G4int A, G4double E,
                       G4Ions::G4FloatLevelBase flb)
{

   G4IsotopeProperty* fProperty = nullptr;

   // At first searching UserDefined
   if ( fUserDefinedList ) {
      for ( G4IsotopeList::iterator it = fUserDefinedList->begin() ; it != fUserDefinedList->end() ; it++ ) {
        
         if ( Z == (*it)->GetAtomicNumber() && A == (*it)->GetAtomicMass() ) {
            G4double levelE = (*it)->GetEnergy();         
            if ( levelE - flevelTolerance/2 <= E && E < levelE + flevelTolerance/2 ) {
               if( flb == (*it)->GetFloatLevelBase() )
               { return *it; } //found 
            }
         }
      }
   } 

   //Serching pre-load
   //Note: isomer level is properly set only for pre_load_list.
   //
   G4int ionCode = 1000*Z + A;
   std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > >::iterator itf = map_pre_load_list.find( ionCode );

   if ( itf !=  map_pre_load_list.end() ) {
      std::multimap< G4double , G4IsotopeProperty* >::iterator lower_bound_itr = itf -> second.lower_bound ( E - flevelTolerance/2 );
      G4double levelE = DBL_MAX;
/*
      if ( lower_bound_itr !=  itf -> second.end() ) {
         levelE = lower_bound_itr->first;
         if ( levelE - flevelTolerance/2 <= E && E < levelE + flevelTolerance/2 ) {
            if( flb == (lower_bound_itr->second)->GetFloatLevelBase() )
            { return lower_bound_itr->second; } // found 
         }
      }
*/
      while ( lower_bound_itr != itf -> second.end() ) {
         levelE = lower_bound_itr->first;
         if ( levelE - flevelTolerance/2 <= E && E < levelE + flevelTolerance/2 ) {
            if ( flb == (lower_bound_itr->second)->GetFloatLevelBase() ) return lower_bound_itr->second; // found 
         } else {
            break;
         } 
         lower_bound_itr++;    
      }
        
   }

   return fProperty; // not found
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
   ;
}

///////////////////////////////////////////////////////////////////////////////
void G4NuclideTable::GenerateNuclide()
{

   if ( threshold_of_half_life < minimum_threshold_of_half_life ) {

      // Need to update full list

      char* path = getenv("G4ENSDFSTATEDATA");

      if ( !path ) {
         G4Exception("G4NuclideTable", "PART70000",
                     FatalException, "G4ENSDFSTATEDATA environment variable must be set");
         return;
      }

      std::ifstream ifs;
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
      G4String ionFL;
      G4double ionLife;
      G4int ionJ;
      G4double ionMu;

      //ifs >> ionZ >> ionA >> ionE >> ionLife >> ionJ >> ionMu;
      ifs >> ionZ >> ionA >> ionE >> ionFL >> ionLife >> ionJ >> ionMu;

      while ( ifs.good() ) {// Loop checking, 09.08.2015, K.Kurashige

         if ( ionCode != 1000*ionZ + ionA ) {
              iLevel = 0;
              ionCode = 1000*ionZ + ionA;
         }

         ionE *= keV;
         //G4int flbIndex = 0;
         //ionE = StripFloatLevelBase( ionE, flbIndex );
         G4Ions::G4FloatLevelBase flb = StripFloatLevelBase( ionFL );
         ionLife *= ns;
         ionMu *= (joule/tesla);

         //if ( ( ionE == 0 && minimum_threshold_of_half_life == DBL_MAX ) // ground state is alwyas build in very first attempt
         //  || ( threshold_of_half_life <= ionLife*std::log(2.0) && ionLife*std::log(2.0) < minimum_threshold_of_half_life && ionE != 0 ) ) {

         if ( ( ionE == 0 && flb == G4Ions::G4FloatLevelBase::no_Float ) //
           || ( threshold_of_half_life <= ionLife*std::log(2.0) && ionLife*std::log(2.0) < minimum_threshold_of_half_life ) ) {

            if ( ionE > 0 ) iLevel++;
            if ( iLevel > 9 ) iLevel=9;

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
            fProperty->SetFloatLevelBase( flb );

            fIsotopeList->push_back(fProperty);

            std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > >::iterator itf = map_full_list.find( ionCode );
            if ( itf == map_full_list.end() ) {
               std::multimap<G4double, G4IsotopeProperty*> aMultiMap;
               //itf = map_full_list.insert( std::pair< G4int , std::multimap< G4double , G4IsotopeProperty* > > ( ionCode , aMultiMap ) );
               itf = ( map_full_list.insert( std::pair< G4int , std::multimap< G4double , G4IsotopeProperty* > > ( ionCode , aMultiMap ) ) ).first;
            } 
            itf -> second.insert( std::pair< G4double, G4IsotopeProperty* >( ionE , fProperty ) );
         }

         ifs >> ionZ >> ionA >> ionE >> ionFL >> ionLife >> ionJ >> ionMu;
      }

      minimum_threshold_of_half_life = threshold_of_half_life;

   }


   // Clear current map
   for ( std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > >::iterator
      it = map_pre_load_list.begin(); it != map_pre_load_list.end(); it++ ) {
      it->second.clear();
   }
   map_pre_load_list.clear();

   // Build map based on current threshold value 
   for ( std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > >::iterator
      it = map_full_list.begin(); it != map_full_list.end(); it++ ) {

      G4int ionCode = it->first;
      std::map< G4int , std::multimap< G4double , G4IsotopeProperty* > >::iterator itf = map_pre_load_list.find( ionCode );
      if ( itf == map_pre_load_list.end() ) {
         std::multimap<G4double, G4IsotopeProperty*> aMultiMap;
         itf = ( map_pre_load_list.insert( std::pair< G4int , std::multimap< G4double , G4IsotopeProperty* > > ( ionCode , aMultiMap ) ) ).first;
      }
      G4int iLevel = 0;
      for ( std::multimap< G4double , G4IsotopeProperty* >::iterator
         itt = it->second.begin(); itt != it->second.end(); itt++ ) {

         G4double exEnergy = itt->first;
         G4double meanLife = itt->second->GetLifeTime();

         if ( exEnergy == 0.0
           || meanLife*std::log(2.0) > threshold_of_half_life ) {

            if ( itt->first != 0.0 ) iLevel++;
            if ( iLevel > 9 ) iLevel=9;
            itt->second->SetIsomerLevel( iLevel );

            itf -> second.insert( std::pair< G4double, G4IsotopeProperty* >( exEnergy , itt->second ) );
         }
      }
   }

}

void G4NuclideTable::AddState( G4int ionZ, G4int ionA, G4double ionE, G4double ionLife, G4int ionJ, G4double ionMu )
{
   if ( G4Threading::IsMasterThread() ) {
      G4int flbIndex = 0;
      ionE = StripFloatLevelBase( ionE, flbIndex );
      AddState(ionZ,ionA,ionE,flbIndex,ionLife,ionJ,ionMu);
   }
}

void G4NuclideTable::AddState( G4int ionZ, G4int ionA, G4double ionE,
                               G4int flbIndex, G4double ionLife, G4int ionJ, G4double ionMu )
{
   if ( G4Threading::IsMasterThread() ) {

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
   fProperty->SetFloatLevelBase(flbIndex);

   fUserDefinedList->push_back(fProperty);
   fIsotopeList->push_back(fProperty);

   }
}

void G4NuclideTable::AddState( G4int ionZ, G4int ionA, G4double ionE,
                               G4Ions::G4FloatLevelBase flb, G4double ionLife, G4int ionJ, G4double ionMu )
{
   if ( G4Threading::IsMasterThread() ) {

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
      fProperty->SetFloatLevelBase( flb );

      fUserDefinedList->push_back(fProperty);
      fIsotopeList->push_back(fProperty);

   }
}

//#include "G4Threading.hh"
void G4NuclideTable::SetThresholdOfHalfLife( G4double t )
{
   if ( G4Threading::IsMasterThread() ) {
      threshold_of_half_life=t; 
      GenerateNuclide();
   }
}

G4double G4NuclideTable::StripFloatLevelBase(G4double E, G4int& flbIndex)
{
  G4double rem = std::fmod(E/(1.0E-3*eV),10.0);
  flbIndex = int(rem);
  return E-rem;
}

G4Ions::G4FloatLevelBase G4NuclideTable::StripFloatLevelBase( G4String sFLB )
{
   if ( sFLB.size() < 1 || 2 < sFLB.size() ) {
      G4String text;
      text += sFLB; 
      text += " is not valid indicator of G4Ions::G4FloatLevelBase. You may use a wrong version of ENSDFSTATE data. Please use G4ENSDFSTATE2.0 or later.";

      G4Exception( "G4NuclideTable" , "PART70002" ,
                   FatalException , text );
   }
   G4Ions::G4FloatLevelBase flb = noFloat;
   if ( !(sFLB == '-') ) {
      flb = G4Ions::FloatLevelBase( sFLB.back() );
   }
   return flb;
}
