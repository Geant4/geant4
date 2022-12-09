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
// G4NuclideTable class implementation
//
// Author: T.Koi, SLAC - 10 October 2013
// --------------------------------------------------------------------

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
 
// --------------------------------------------------------------------
G4NuclideTable* G4NuclideTable::GetInstance()
{
  static G4NuclideTable instance;
  return &instance;
}

// --------------------------------------------------------------------
G4NuclideTable* G4NuclideTable::GetNuclideTable()
{
  return GetInstance();
}

// --------------------------------------------------------------------
G4NuclideTable::G4NuclideTable()
  : G4VIsotopeTable("Isomer"),
    mean_life_threshold(1.0*ns),
    flevelTolerance(1.0*eV)
{
  fMessenger = new G4NuclideTableMessenger(this);
  fIsotopeList = new G4IsotopeList();
  GenerateNuclide();
}

// --------------------------------------------------------------------
G4NuclideTable::~G4NuclideTable()
{
  for (auto it=map_pre_load_list.begin(); it!=map_pre_load_list.end(); ++it)
  {
    it->second.clear();
  }
  map_pre_load_list.clear();

  for (auto it=map_full_list.begin(); it!=map_full_list.end(); ++it)
  {
    it->second.clear();
  }
  map_full_list.clear();

  if (fIsotopeList != nullptr)
  {
    for (std::size_t i = 0 ; i<fIsotopeList->size(); ++i)
    {
      delete (*fIsotopeList)[i];
    }
    fIsotopeList->clear();
    delete fIsotopeList;
    fIsotopeList = nullptr;
  }
  delete fMessenger;
}

// --------------------------------------------------------------------
G4IsotopeProperty* G4NuclideTable::GetIsotope(G4int Z, G4int A, G4double E,
                                              G4Ions::G4FloatLevelBase flb)
{
  G4IsotopeProperty* fProperty = nullptr;

  // At first searching UserDefined
  if ( fUserDefinedList )
  {
    for (auto it=fUserDefinedList->cbegin(); it!=fUserDefinedList->cend(); ++it)
    {
      if ( Z == (*it)->GetAtomicNumber() && A == (*it)->GetAtomicMass() )
      {
        G4double levelE = (*it)->GetEnergy();         
        if ( levelE - flevelTolerance/2 <= E && E < levelE + flevelTolerance/2 )
        {
          if( flb == (*it)->GetFloatLevelBase() ) { return *it; } //found 
        }
      }
    }
  } 

  // Searching pre-load
  // Note: isomer level is properly set only for pre_load_list
  //
  G4int ionCode = 1000*Z + A;
  auto itf = map_pre_load_list.find( ionCode );

  if ( itf !=  map_pre_load_list.cend() )
  {
    auto lower_bound_itr = itf -> second.lower_bound ( E - flevelTolerance/2 );
    G4double levelE = DBL_MAX;

    while ( lower_bound_itr != itf -> second.cend() )
    {
      levelE = lower_bound_itr->first;
      if ( levelE - flevelTolerance/2 <= E && E < levelE + flevelTolerance/2 )
      {
        if ( flb == (lower_bound_itr->second)->GetFloatLevelBase() )
        {
          return lower_bound_itr->second; // found
        }
      }
      else
      {
        break;
      } 
      ++lower_bound_itr;    
    }
  }

  return fProperty; // not found
}

// --------------------------------------------------------------------
G4double G4NuclideTable::GetTruncationError( G4double eex )
{
  G4double tolerance= G4NuclideTable::GetInstance()->GetLevelTolerance();
  return eex - (G4long)(eex/tolerance)*tolerance;
}

// --------------------------------------------------------------------
G4double G4NuclideTable::Round( G4double eex )
{
  G4double tolerance= G4NuclideTable::GetInstance()->GetLevelTolerance();
  return round(eex/tolerance)*tolerance;
}

// --------------------------------------------------------------------
G4long G4NuclideTable::Truncate( G4double eex )
{
  G4double tolerance= G4NuclideTable::GetInstance()->GetLevelTolerance();
  return (G4long)(eex/tolerance);
}

// --------------------------------------------------------------------
G4double G4NuclideTable::Tolerance()
{
  return G4NuclideTable::GetInstance()->GetLevelTolerance();
}

// --------------------------------------------------------------------
G4IsotopeProperty*
G4NuclideTable::GetIsotopeByIsoLvl(G4int Z, G4int A, G4int lvl)
{
  if(lvl==0) return GetIsotope(Z,A,0.0);
  return nullptr;
}

// --------------------------------------------------------------------
void G4NuclideTable::GenerateNuclide()
{
  if (mean_life_threshold < minimum_mean_life_threshold) {
    // Need to update full list
    const char* path = G4FindDataDir("G4ENSDFSTATEDATA");

    if (path == nullptr) {
      G4Exception("G4NuclideTable", "PART70000", FatalException,
                  "G4ENSDFSTATEDATA environment variable must be set");
      return;
    }

    std::ifstream ifs;
    G4String filename(path);
    filename += "/ENSDFSTATE.dat";

    ifs.open(filename.c_str() );
    if (!ifs.good() ) {
      G4Exception("G4NuclideTable", "PART70001", FatalException,
                  "ENSDFSTATE.dat is not found.");
      return;
    }

    G4int ionCode = 0;
    G4int iLevel = 0;
    G4int ionZ;
    G4int ionA;
    G4double ionE;
    G4String ionFL;
    G4double ionLife;
    G4int ionJ;
    G4double ionMu;

    // Lifetimes read from ENSDFSTATE are mean lives
    ifs >> ionZ >> ionA >> ionE >> ionFL >> ionLife >> ionJ >> ionMu;

    while (ifs.good() )  // Loop checking, 09.08.2015, K.Kurashige
    {
      if (ionCode != 1000*ionZ + ionA ) {
        iLevel = 0;
        ionCode = 1000*ionZ + ionA;
      }

      ionE *= keV;
      G4Ions::G4FloatLevelBase flb = StripFloatLevelBase(ionFL);
      ionLife *= ns;
      ionMu *= (joule/tesla);

      if ((ionE == 0 && flb == G4Ions::G4FloatLevelBase::no_Float) ||
          (mean_life_threshold <= ionLife &&
           ionLife < minimum_mean_life_threshold) ) {
        if (ionE > 0) ++iLevel;
        if (iLevel > 9) iLevel = 9;

        G4IsotopeProperty* fProperty = new G4IsotopeProperty();

        // Set Isotope Property
        fProperty->SetAtomicNumber(ionZ);
        fProperty->SetAtomicMass(ionA);
        fProperty->SetIsomerLevel(iLevel);
        fProperty->SetEnergy(ionE);
        fProperty->SetiSpin(ionJ);
        fProperty->SetLifeTime(ionLife);
        fProperty->SetDecayTable(nullptr);
        fProperty->SetMagneticMoment(ionMu);
        fProperty->SetFloatLevelBase( flb );

        fIsotopeList->push_back(fProperty);

        auto itf = map_full_list.find(ionCode);
        if (itf == map_full_list.cend() ) {
          std::multimap<G4double, G4IsotopeProperty*> aMultiMap;
          itf = (map_full_list.insert(
                 std::pair<G4int, std::multimap<G4double,
                                  G4IsotopeProperty*> > (ionCode, aMultiMap) ) ).first;
        }
        itf->second.insert(
             std::pair<G4double, G4IsotopeProperty*>(ionE, fProperty) );
      }

      ifs >> ionZ >> ionA >> ionE >> ionFL >> ionLife >> ionJ >> ionMu;
    }  // End while

    minimum_mean_life_threshold = mean_life_threshold;
  }

  // Clear current map
  for (auto it = map_pre_load_list.begin(); it != map_pre_load_list.end(); ++it) {
    it->second.clear();
  }
  map_pre_load_list.clear();

  // Build map based on current threshold value
  for (auto it = map_full_list.cbegin(); it != map_full_list.cend(); ++it) {
    G4int ionCode = it->first;
    auto itf = map_pre_load_list.find(ionCode);
    if (itf == map_pre_load_list.cend() ) {
      std::multimap<G4double, G4IsotopeProperty*> aMultiMap;
      itf = (map_pre_load_list.insert(
             std::pair<G4int, std::multimap<G4double,
                              G4IsotopeProperty*> > (ionCode, aMultiMap) ) ).first;
    }

    G4int iLevel = 0;
    for (auto itt = it->second.cbegin(); itt != it->second.cend(); ++itt) {
      G4double exEnergy = itt->first;
      G4double meanLife = itt->second->GetLifeTime();
      if (exEnergy == 0.0 || meanLife > mean_life_threshold) {
        if (itt->first != 0.0) ++iLevel;
        if (iLevel > 9) iLevel = 9;
        itt->second->SetIsomerLevel(iLevel);
        itf->second.insert(
             std::pair<G4double, G4IsotopeProperty*>(exEnergy, itt->second) );
      }
    }
  }
}

// --------------------------------------------------------------------
void G4NuclideTable::AddState( G4int ionZ, G4int ionA, G4double ionE,
                               G4double ionLife, G4int ionJ, G4double ionMu )
{
  if ( G4Threading::IsMasterThread() )
  {
    G4int flbIndex = 0;
    ionE = StripFloatLevelBase( ionE, flbIndex );
    AddState(ionZ,ionA,ionE,flbIndex,ionLife,ionJ,ionMu);
  }
}

// --------------------------------------------------------------------
void G4NuclideTable::AddState( G4int ionZ, G4int ionA, G4double ionE,
                               G4int flbIndex, G4double ionLife, G4int ionJ,
                               G4double ionMu )
{
  if ( G4Threading::IsMasterThread() )
  {
    if ( fUserDefinedList == nullptr ) fUserDefinedList = new G4IsotopeList();

    G4IsotopeProperty* fProperty = new G4IsotopeProperty(); 

    // Set Isotope Property
    fProperty->SetAtomicNumber(ionZ);
    fProperty->SetAtomicMass(ionA);
    fProperty->SetIsomerLevel(9);
    fProperty->SetEnergy(ionE);
    fProperty->SetiSpin(ionJ);
    fProperty->SetLifeTime(ionLife);
    fProperty->SetDecayTable(nullptr);
    fProperty->SetMagneticMoment(ionMu);
    fProperty->SetFloatLevelBase(flbIndex);

    fUserDefinedList->push_back(fProperty);
    fIsotopeList->push_back(fProperty);
   }
}

// --------------------------------------------------------------------
void G4NuclideTable::AddState( G4int ionZ, G4int ionA, G4double ionE,
                               G4Ions::G4FloatLevelBase flb, G4double ionLife,
                               G4int ionJ, G4double ionMu )
{
  if ( G4Threading::IsMasterThread() )
  {
    if ( fUserDefinedList == nullptr ) fUserDefinedList = new G4IsotopeList();

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
    fProperty->SetFloatLevelBase(flb);

    fUserDefinedList->push_back(fProperty);
    fIsotopeList->push_back(fProperty);
  }
}

// --------------------------------------------------------------------
void G4NuclideTable::SetThresholdOfHalfLife(G4double t) {
  if (G4Threading::IsMasterThread() ) {
    mean_life_threshold = t/0.69314718;
    GenerateNuclide();
  }
}

// Set the mean life threshold for nuclides
// All nuclides with mean lives greater than this value are created 
// for this run  
void G4NuclideTable::SetMeanLifeThreshold(G4double t)
{
  if (G4Threading::IsMasterThread() )
  {
    mean_life_threshold = t;
    GenerateNuclide();
  }
}

// --------------------------------------------------------------------
G4double G4NuclideTable::StripFloatLevelBase(G4double E, G4int& flbIndex)
{
  G4double rem = std::fmod(E/(1.0E-3*eV),10.0);
  flbIndex = G4int(rem);
  return E-rem;
}

// --------------------------------------------------------------------
G4Ions::G4FloatLevelBase
G4NuclideTable::StripFloatLevelBase( const G4String& sFLB )
{
   if ( sFLB.size() < 1 || 2 < sFLB.size() )
   {
     G4String text;
     text += sFLB; 
     text += " is not valid indicator of G4Ions::G4FloatLevelBase.\n";
     text += "You may use a wrong version of ENSDFSTATE data.\n";
     text += "Please use G4ENSDFSTATE-2.0 or later.";

     G4Exception( "G4NuclideTable", "PART70002", FatalException, text );
   }
   G4Ions::G4FloatLevelBase flb = noFloat;
   if ( !(sFLB == "-") )
   {
     flb = G4Ions::FloatLevelBase( sFLB.back() );
   }
   return flb;
}
