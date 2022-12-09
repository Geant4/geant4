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
// G4IonTable class implementation
//
// Author: H.Kurashige, 27 June 1998
// --------------------------------------------------------------------

#include <iostream>               
#include <iomanip>               
#include <sstream>
#include <algorithm>
#include <vector>     

#include "G4ios.hh"
#include "G4Threading.hh"
#include "G4AutoDelete.hh"

#include "G4IonTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4StateManager.hh"
#include "G4Ions.hh"
#include "G4UImanager.hh"
#include "G4NucleiProperties.hh"
#include "G4HyperNucleiProperties.hh"

#include "G4IsotopeProperty.hh"
#include "G4VIsotopeTable.hh"
#include "G4NuclideTable.hh"

#include "G4MuonicAtom.hh"
#include "G4MuonicAtomHelper.hh"

// It is very important for multithreaded Geant4 to keep only one copy of the
// particle table pointer and the ion table pointer. However, we try to let 
// each worker thread hold its own copy of the particle dictionary and the 
// ion list. This implementation is equivalent to make the ion table thread
// private. The two shadow ponters are used by each worker thread to copy the
// content from the master thread.
//
G4ThreadLocal G4IonTable::G4IonList* G4IonTable::fIonList = nullptr;
G4ThreadLocal std::vector<G4VIsotopeTable*> *G4IonTable::fIsotopeTableList = nullptr;
G4IonTable::G4IonList* G4IonTable::fIonListShadow = nullptr;
std::vector<G4VIsotopeTable*> *G4IonTable::fIsotopeTableListShadow = nullptr;

namespace lightions
{
  static const G4ParticleDefinition* p_proton = nullptr;
  static const G4ParticleDefinition* p_deuteron = nullptr;
  static const G4ParticleDefinition* p_triton = nullptr;
  static const G4ParticleDefinition* p_alpha = nullptr;
  static const G4ParticleDefinition* p_He3 = nullptr;
  void Init()
  {
    if ( p_proton != nullptr ) return;
    p_proton   = G4ParticleTable::GetParticleTable()->FindParticle("proton");
    p_deuteron = G4ParticleTable::GetParticleTable()->FindParticle("deuteron");
    p_triton   = G4ParticleTable::GetParticleTable()->FindParticle("triton");
    p_alpha    = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    p_He3      = G4ParticleTable::GetParticleTable()->FindParticle("He3");
  }
}

namespace antilightions
{
  static const G4ParticleDefinition* p_proton = nullptr;
  static const G4ParticleDefinition* p_deuteron = nullptr;
  static const G4ParticleDefinition* p_triton = nullptr;
  static const G4ParticleDefinition* p_alpha = nullptr;
  static const G4ParticleDefinition* p_He3 = nullptr;
  void Init()
  {
    if ( p_proton != nullptr ) return;
    p_proton   = G4ParticleTable::GetParticleTable()->FindParticle("anti_proton");
    p_deuteron = G4ParticleTable::GetParticleTable()->FindParticle("anti_deuteron");
    p_triton   = G4ParticleTable::GetParticleTable()->FindParticle("anti_triton");
    p_alpha    = G4ParticleTable::GetParticleTable()->FindParticle("anti_alpha");
    p_He3      = G4ParticleTable::GetParticleTable()->FindParticle("anti_He3");
  }
}

#ifdef G4MULTITHREADED
G4Mutex G4IonTable::ionTableMutex = G4MUTEX_INITIALIZER;
#endif 

// --------------------------------------------------------------------
// Constructor
//
G4IonTable::G4IonTable()
{
  fIonList = new G4IonList();

  // Set up the shadow pointer used by worker threads.
  //
  if (fIonListShadow == nullptr)
  {
    fIonListShadow = fIonList;
  }

  fIsotopeTableList = new std::vector<G4VIsotopeTable*>;

  // Set up the shadow pointer used by worker threads.
  //
  if (fIsotopeTableListShadow == nullptr)
  {
    fIsotopeTableListShadow = fIsotopeTableList;
  }    

  PrepareNuclideTable();
  RegisterIsotopeTable(pNuclideTable);
}

// --------------------------------------------------------------------
// Destructor
//
G4IonTable::~G4IonTable()
{
  // delete IsotopeTable if existing
  if (fIsotopeTableList !=  nullptr )
  {
    for (std::size_t i=0; i<fIsotopeTableList->size(); ++i)
    {
      G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[i];
      if( fIsotopeTable != G4NuclideTable::GetNuclideTable() )
      {
        delete fIsotopeTable;
      }
    }
    fIsotopeTableList->clear();
    delete fIsotopeTableList;
  }
  fIsotopeTableList = nullptr;

  if (fIonList == nullptr) return;

  // remove all contents in the Ion List 
  // No need to delete here because all particles are dynamic objects
  fIonList->clear();
  delete fIonList;
  fIonList = nullptr;
}

// --------------------------------------------------------------------
// GetIonTable
//
G4IonTable* G4IonTable::GetIonTable()
{
  return G4ParticleTable::GetParticleTable()->GetIonTable();
}

// --------------------------------------------------------------------
// WorkerG4IonTable
//
// Used by each worker thread to copy the content from the master thread
//
void G4IonTable::WorkerG4IonTable()
{
  if( fIonList == nullptr ) { fIonList = new G4IonList(); }
  else                      { fIonList->clear(); }

  for (auto it = fIonListShadow->cbegin(); it != fIonListShadow->cend(); ++it )
  {
    fIonList->insert(*it);
  }

  // Do not copy Isotope Table to Worker thread
  //
  if( fIsotopeTableList == nullptr )
  {
    fIsotopeTableList = new std::vector<G4VIsotopeTable*>; 
    for (std::size_t i = 0; i < fIsotopeTableListShadow->size(); ++i)
    {
      fIsotopeTableList->push_back((*fIsotopeTableListShadow)[i]); 
    }
  }
}

// --------------------------------------------------------------------
// InitializeLightIons
//
void G4IonTable::InitializeLightIons()
{
  lightions::Init();
  antilightions::Init();
}

// --------------------------------------------------------------------
// DestroyWorkerG4IonTable
//
void G4IonTable::DestroyWorkerG4IonTable()
{
  // delete IsotopeTable if existing
  if (fIsotopeTableList != nullptr )
  {
    for (std::size_t i=0; i<fIsotopeTableList->size(); ++i)
    {
      G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[i];
      if( fIsotopeTable != G4NuclideTable::GetNuclideTable() )
      {
        delete fIsotopeTable;
      }
    }
    fIsotopeTableList->clear();
    delete fIsotopeTableList;
  }
  fIsotopeTableList = nullptr;


  if (fIonList == nullptr) return;

  // remove all contents in the Ion List 
  // No need to delete here because all particles are dynamic objects
  fIonList->clear();
  delete fIonList;
  fIonList = nullptr;
}

// --------------------------------------------------------------------
// CreateIon
//
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, G4double E, 
                                            G4Ions::G4FloatLevelBase flb)
{
  G4ParticleDefinition* ion = nullptr;

  // check whether GenericIon has processes
  G4ParticleDefinition* genericIon = 
    G4ParticleTable::GetParticleTable()->GetGenericIon();
  G4ProcessManager* pman = nullptr;
  if (genericIon!= nullptr) { pman = genericIon->GetProcessManager(); }
  if ((genericIon == nullptr) || (genericIon->GetParticleDefinitionID() < 0)
                              || (pman==nullptr))
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1)
    {
      G4cout << "G4IonTable::CreateIon() : can not create ion of  " 
             << " Z =" << Z << "  A = " << A 
             << "  because GenericIon is not ready !!" << G4endl;
    }
#endif
    G4Exception( "G4IonTable::CreateIon()", "PART105", JustWarning, 
                 "Can not create ions because GenericIon is not ready");
    return nullptr;
  }
  
  G4double life = 0.0;
  G4DecayTable* decayTable = nullptr;
  G4bool stable = true;
  G4double mu  = 0.0;
  G4double Eex = 0.0;
  G4int    lvl = 0;
  G4int    J   = 0;

  const G4IsotopeProperty* fProperty = FindIsotope(Z, A, E, flb);
  if (fProperty != nullptr )
  {
    Eex  = fProperty->GetEnergy();
    lvl  = fProperty->GetIsomerLevel();
    J    = fProperty->GetiSpin();
    life = fProperty->GetLifeTime();
    mu   = fProperty->GetMagneticMoment();    
    decayTable = fProperty->GetDecayTable();
    stable = (life <= 0.) || (decayTable == nullptr);
    lvl = fProperty->GetIsomerLevel();
    if (lvl <0) lvl=9;
  }
  else
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1)
    {
      G4ExceptionDescription ed;
      ed << "G4IonTable::CreateIon(): G4IsotopeProperty object is not found for"
         << " Z = " << Z << " A = " << A << " E = " << E/keV << " (keV)";
      if(flb!=G4Ions::G4FloatLevelBase::no_Float)
      {
        ed << " FloatingLevel +" << G4Ions::FloatLevelBaseChar(flb); 
      }
      ed << ".\n"
         << " Physics quantities such as life are not set for this ion.";
      G4Exception( "G4IonTable::CreateIon()","PART70105", JustWarning, ed);
    }
#endif
    // excitation energy
    Eex = E;
    // lvl is assigned to 9 temporarily    
    if (Eex>0.0) lvl=9;
  }

  // Eex = G4NuclideTable::Round(Eex); 
  if (Eex==0.0) lvl=0;
  // ion name
  G4String name =""; 
  /////////////if (lvl<9) name = GetIonName(Z, A, lvl);
  if (lvl==0 && flb==G4Ions::G4FloatLevelBase::no_Float)
    name = GetIonName(Z, A, lvl);
  else
    name = GetIonName(Z, A, Eex, flb);

  // PDG encoding 
  G4int encoding = GetNucleusEncoding(Z,A,E,lvl);

  // PDG mass
  G4double mass =  GetNucleusMass(Z, A)+ Eex;
 
  // PDG charge is set to one of nucleus
  G4double charge =  G4double(Z)*eplus;
 
  // create an ion
  // spin, parity, isospin values are fixed

  // Request lock for particle table accesses. Some changes are inside 
  // this critical region.
  //
  ion = new G4Ions(   name,            mass,       0.0*MeV,     charge, 
                         J,              +1,             0,          
                         0,               0,             0,             
                 "nucleus",               0,             A,    encoding,
                    stable,            life,    decayTable,       false,
                 "generic",               0,
                       Eex,             lvl         );

  // Release lock for particle table accesses.
  //
  ion->SetPDGMagneticMoment(mu);
  static_cast<G4Ions*>(ion)->SetFloatLevelBase(flb);

  // No Anti particle registered
  ion->SetAntiPDGEncoding(0);
  
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
  {
    G4cout << "G4IonTable::CreateIon() : create ion of " << name
           << "  " << Z << ", " << A
           << " encoding=" << encoding;
    if (E>0.0)
    {
      G4cout << " IsomerLVL=" << lvl
             << " excited energy=" << Eex/keV << "[keV]";
    }
    G4cout << G4endl;
  } 
#endif
  
  // Add process manager to the ion
  AddProcessManager(ion);
 
#ifdef G4MULTITHREADED
  // Fill decay channels if this method is invoked from worker
  if(G4Threading::IsWorkerThread())
  {
    if(!stable && decayTable)
    {
      G4int nCh = decayTable->entries();
      for(G4int iCh=0; iCh<nCh; ++iCh)
      { 
        decayTable->GetDecayChannel(iCh)->GetDaughter(0); 
      }
    }
  }
#endif
  
  return ion;
}

// --------------------------------------------------------------------
// CreateIon
//
G4ParticleDefinition*
G4IonTable::CreateIon(G4int Z, G4int A, G4int LL, G4double E,
                      G4Ions::G4FloatLevelBase flb)
{
  if (LL==0) return CreateIon(Z,A,E,flb);
  
  // create hyper nucleus
  G4ParticleDefinition* ion = nullptr;

  // check whether GenericIon has processes
  G4ParticleDefinition* genericIon = 
    G4ParticleTable::GetParticleTable()->GetGenericIon();
  G4ProcessManager* pman = nullptr;
  if (genericIon != nullptr) pman = genericIon->GetProcessManager();
  if ((genericIon == nullptr) || (genericIon->GetParticleDefinitionID() < 0)
                              || (pman==nullptr))
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1)
    {
      G4cout << "G4IonTable::CreateIon() : can not create ion of  " 
             << " Z =" << Z << "  A = " << A 
             << "  because GenericIon is not ready !!" << G4endl;
    }
#endif
    G4Exception( "G4IonTable::CreateIon()","PART105", JustWarning, 
                 "Can not create ions because GenericIon is not ready");
    return nullptr;
  }
 
  G4int J = 0;
  G4double life = 0.0;
  G4DecayTable* decayTable = nullptr;
  G4bool stable = true;
 
  // excitation energy
  // G4double Eex = G4NuclideTable::Round(E); 
  G4double Eex = E; 
  G4double mass =  GetNucleusMass(Z, A, LL)+ Eex;
  G4int    lvl = 0;
  // lvl is assigned to 9 temporarily
  if (Eex>0.0) lvl=9;
   
  // PDG encoding 
  G4int encoding = GetNucleusEncoding(Z,A,LL,E,lvl);

  // PDG charge is set to one of nucleus
  G4double charge =  G4double(Z)*eplus;

  // create an ion
  // spin, parity, isospin values are fixed
  //
  // get ion name
  G4String name = GetIonName(Z, A, LL, Eex, flb);
  
  ion = new G4Ions(   name,            mass,       0.0*MeV,     charge, 
                         J,              +1,             0,          
                         0,               0,             0,             
                 "nucleus",               0,             A,    encoding,
                    stable,            life,    decayTable,       false,
                  "generic",              0,
                      Eex,              lvl         );

  // Release lock for particle table accesses
 
  G4double mu = 0.0; //  magnetic moment
  ion->SetPDGMagneticMoment(mu);
  static_cast<G4Ions*>(ion)->SetFloatLevelBase(flb);

  // No Anti particle registered
  ion->SetAntiPDGEncoding(0);
  
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
  {
    G4cout << "G4IonTable::CreateIon() : create hyper ion of " << name
           << "  " << Z << ", " << A << ", " << LL
           << " encoding=" << encoding;
    if (E>0.0)
    {
      G4cout << " IsomerLVL=" << lvl
             << " excited energy=" << Eex/keV << "[keV]";
    }
    G4cout << G4endl;
  } 
#endif
  
  // Add process manager to the ion
  AddProcessManager(ion);
      
  return ion;
}

// --------------------------------------------------------------------
// CreateIon
//
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, G4int lvl)
{
  if(lvl == 0) return CreateIon(Z,A,0.0,G4Ions::G4FloatLevelBase::no_Float);
  G4Exception( "G4IonTable::CreateIon()","PART105", JustWarning, 
      "Ion cannot be created by an isomer level. Use excitation energy.");
  return nullptr;
}

// --------------------------------------------------------------------
// CreateIon
//
G4ParticleDefinition*
G4IonTable::CreateIon(G4int Z, G4int A, G4int LL, G4int lvl)
{
  if (LL==0) return CreateIon(Z,A,lvl);
  if(lvl == 0) return CreateIon(Z,A,0.0,G4Ions::G4FloatLevelBase::no_Float);
  
  if (lvl>0)
  {
    G4ExceptionDescription ed;
    ed << "Isomer level " << lvl << " is unknown for the isotope (Z="
       << Z << ", A=" << A << ", L=" << LL << "). Null pointer is returned.";
    G4Exception( "G4IonTable::GetIon()","PART106", JustWarning, ed);
    return nullptr;
  }
  return nullptr;
}

// --------------------------------------------------------------------
// -- GetIon methods  ------
//
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int lvl)
{
  if ( (A<1) || (Z<=0) || (lvl<0) || (A>999) )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A <<  "  Lvl = " << lvl << G4endl;
    }
#endif
    return nullptr;
  }
  if ( lvl == 0 ) return GetIon(Z,A,0.0);
  
  // Search ions with A, Z, lvl 
  G4ParticleDefinition* ion = FindIon(Z,A,lvl);
  
  // create ion
#ifdef G4MULTITHREADED
  if (ion == nullptr )
  {
    if(G4Threading::IsWorkerThread())
    {
      G4MUTEXLOCK(&G4IonTable::ionTableMutex);
      ion = FindIonInMaster(Z,A,lvl);
      if(ion != nullptr ) InsertWorker(ion); 
      G4MUTEXUNLOCK(&G4IonTable::ionTableMutex);
    } 
  }
#endif
  if (ion == nullptr )
  {
    G4Exception( "G4IonTable::GetIon()","PART105", JustWarning, 
      "Ion cannot be created by an isomer level. Use excitation energy.");
  }
  return ion;  
}

// --------------------------------------------------------------------
//
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int LL, G4int lvl)
{
  if (LL==0) return GetIon(Z,A,lvl);

  if (A < 2 || Z < 0 || Z > A-LL || LL>A || A>999 )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A << " L = " << LL 
             <<"  IsomerLvl = " << lvl << G4endl;
    }
#endif
    return nullptr;
  }
  else if( A==2 )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::GetIon() : No boud state for " 
             << " Z =" << Z << "  A = " << A << " L = " << LL 
             <<"  IsomerLvl = " << lvl << G4endl;
    }
#endif
    return nullptr;
  }

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,LL,lvl);

  // create ion
  if (ion == nullptr)
  {
    if (lvl==0)
    {
#ifdef G4MULTITHREADED
      if(G4Threading::IsWorkerThread())
      {
        G4MUTEXLOCK(&G4IonTable::ionTableMutex);
        ion = FindIonInMaster(Z,A,LL,lvl);
        if(ion == nullptr) ion = CreateIon(Z, A, LL, lvl);
        InsertWorker(ion);
        G4MUTEXUNLOCK(&G4IonTable::ionTableMutex);
      }
      else
      { 
        ion = CreateIon(Z, A, LL, lvl); 
      }
#else
      ion = CreateIon(Z, A, LL, lvl);
#endif
    } 
  }

  return ion;  
}

// --------------------------------------------------------------------
//
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4double E, G4int J)
{ 
  return GetIon(Z,A,E,G4Ions::G4FloatLevelBase::no_Float,J); 
}

// --------------------------------------------------------------------
//
G4ParticleDefinition*
G4IonTable::GetIon(G4int Z, G4int A, G4double E, char flbChar, G4int J)
{ 
  return GetIon(Z,A,E,G4Ions::FloatLevelBase(flbChar),J); 
}

// --------------------------------------------------------------------
//
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4double E,
                          G4Ions::G4FloatLevelBase flb, G4int J)
{
  if ( (A<1) || (Z<=0) || (E<0.0) || (A>999) || (J<0) )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A <<  "  E = " << E/keV << G4endl;
    }
#endif
    return nullptr;
  }

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,E,flb,J);

  // create ion
#ifdef G4MULTITHREADED
  if(ion == nullptr )
  {
    if(G4Threading::IsWorkerThread())
    {
      G4MUTEXLOCK(&G4IonTable::ionTableMutex);
      ion = FindIonInMaster(Z,A,E,flb,J);
      if(ion == nullptr) ion = CreateIon(Z,A,E,flb);
      InsertWorker(ion);
      G4MUTEXUNLOCK(&G4IonTable::ionTableMutex);
    }
    else
    { 
      ion = CreateIon(Z,A,E,flb); 
    }
  }
#else
  if (ion == nullptr) ion = CreateIon(Z,A,E,flb);
#endif

  return ion;  
}

// --------------------------------------------------------------------
//
G4ParticleDefinition*
G4IonTable::GetIon(G4int Z, G4int A, G4int LL, G4double E, G4int J)
{ 
  return GetIon(Z,A,LL,E,G4Ions::G4FloatLevelBase::no_Float,J); 
}

// --------------------------------------------------------------------
//
G4ParticleDefinition*
G4IonTable::GetIon(G4int Z, G4int A, G4int LL, G4double E,
                   char flbChar, G4int J)
{ 
  return GetIon(Z,A,LL,E,G4Ions::FloatLevelBase(flbChar),J); 
}

// --------------------------------------------------------------------
//
G4ParticleDefinition*
G4IonTable::GetIon(G4int Z, G4int A, G4int LL, G4double E,
                   G4Ions::G4FloatLevelBase flb, G4int J)
{
  if (LL==0) return GetIon(Z,A,E,flb,J);

  if (A < 2 || Z < 0 || Z > A-LL || LL>A || A>999 )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A << " L = " << LL 
             <<"  E = " << E/keV << G4endl;
    }
#endif
    return nullptr;
  }
  else if( A==2 )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::GetIon() : No boud state for " 
             << " Z =" << Z << "  A = " << A << " L = " << LL 
             <<  "  E = " << E/keV << G4endl;
    }
#endif
    return nullptr;
  }

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,LL,E,flb,J);

  // create ion
#ifdef G4MULTITHREADED
  if(ion == nullptr ){
    if(G4Threading::IsWorkerThread())
    {
      G4MUTEXLOCK(&G4IonTable::ionTableMutex);
      ion = FindIonInMaster(Z,A,LL,E,flb,J);
      if(ion == nullptr) ion = CreateIon(Z,A,LL,E,flb);
      InsertWorker(ion);
      G4MUTEXUNLOCK(&G4IonTable::ionTableMutex);
    }
    else
    { 
      ion = CreateIon(Z,A,LL,E,flb); 
    }
  }
#else
  if(ion == nullptr) ion = CreateIon(Z,A,LL,E,flb);
#endif

  return ion;  
}

// --------------------------------------------------------------------
//
G4ParticleDefinition* G4IonTable::GetIon(G4int encoding)
{
  G4int Z, A, LL, IsoLvl;
  G4double E;
  if (!GetNucleusByEncoding(encoding,Z,A,LL,E,IsoLvl))
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::GetIon() : illegal encoding" 
             << " CODE:" << encoding << G4endl;
    }
#endif
    G4Exception( "G4IonTable::GetIon()","PART106",
                 JustWarning, "illegal encoding for an ion");
    return nullptr;
  }
  return GetIon( Z, A, LL, IsoLvl);
}

// --------------------------------------------------------------------
// -- FindIon methods  ------
//
G4ParticleDefinition*
G4IonTable::FindIon(G4int Z, G4int A, G4double E, G4int J)
{ 
  return FindIon(Z,A,E,G4Ions::G4FloatLevelBase::no_Float,J); 
}

// --------------------------------------------------------------------
//
G4ParticleDefinition*
G4IonTable::FindIon(G4int Z, G4int A, G4double E, char flbChar, G4int J)
{ 
  return FindIon(Z,A,E,G4Ions::FloatLevelBase(flbChar),J); 
}

// --------------------------------------------------------------------
//
G4ParticleDefinition*
G4IonTable::FindIon(G4int Z, G4int A, G4double E,
                    G4Ions::G4FloatLevelBase flb, G4int J)
{
  if ( (A<1) || (Z<=0) || (J<0) || (E<0.0) || (A>999) )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::FindIon(): illegal atomic number/mass"
             << " or excitation level:" << G4endl
             << " Z =" << Z << "  A = " << A <<  "  E = " << E/keV << G4endl;
    }
#endif
    G4Exception("G4IonTable::FindIon()","PART107",
                JustWarning, "illegal atomic number/mass");
    return nullptr;
  }
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion = nullptr;
  G4bool isFound = false;

  // check if light ion
  ion = GetLightIon(Z,A);
  if (ion!= nullptr && E == 0.0)
  { 
    // light ion 
    isFound = true;
  }
  else
  {    
    // -- loop over all particles in Ion table
    G4int encoding = GetNucleusEncoding(Z, A);
    for(auto i = fIonList->find(encoding); i != fIonList->cend(); ++i)
    {
      ion = i->second;
      if ( (ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
      // excitation level
      G4double anExcitaionEnergy= ((const G4Ions*)(ion))->GetExcitationEnergy();
      if (std::fabs(E - anExcitaionEnergy) < pNuclideTable->GetLevelTolerance())
      {
        if(((const G4Ions*)(ion))->GetFloatLevelBase()==flb)
        {
          isFound = true;
          break;
        }
      }
    }
  }

  if ( isFound )
  { 
    return const_cast<G4ParticleDefinition*>(ion);
  }
  else
  {
    return nullptr;
  }
}

// --------------------------------------------------------------------
//
G4ParticleDefinition*
G4IonTable::FindIon(G4int Z, G4int A, G4int LL, G4double E, G4int J)
{ 
  return FindIon(Z,A,LL,E,G4Ions::G4FloatLevelBase::no_Float,J); 
}

// --------------------------------------------------------------------
//
G4ParticleDefinition*
G4IonTable::FindIon(G4int Z, G4int A, G4int LL, G4double E,
                    char flbChar, G4int J)
{ 
  return FindIon(Z,A,LL,E,G4Ions::FloatLevelBase(flbChar),J); 
}

// --------------------------------------------------------------------
//
G4ParticleDefinition*
G4IonTable::FindIon(G4int Z, G4int A, G4int LL, G4double E,
                    G4Ions::G4FloatLevelBase flb, G4int J)
{
  if (LL==0) return FindIon(Z,A,E,flb,J);
  
  if (A < 2 || Z < 0 || Z > A-LL || LL>A || A>999 )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::FindIon(): illegal atomic number/mass"
             << " or excitation level:" << G4endl
             << " Z =" << Z << "  A = " << A << " L = " << LL 
             <<"  E = " << E/keV << G4endl;
    }    
#endif
    G4Exception("G4IonTable::FindIon()", "PART107",
                JustWarning, "illegal atomic number/mass");
    return nullptr;
  }
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion = nullptr;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A, LL, 0.0, 0);
  for(auto i = fIonList->find(encoding); i != fIonList->cend() ; ++i)
  {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if( ion->GetQuarkContent(3) != LL ) break;
    // excitation level
    G4double anExcitaionEnergy = ((const G4Ions*)(ion))->GetExcitationEnergy();
    if (std::fabs(E - anExcitaionEnergy) < pNuclideTable->GetLevelTolerance())
    {
      if(((const G4Ions*)(ion))->GetFloatLevelBase()==flb)
      {
        isFound = true;
        break;
      }
    }
  }

  if ( isFound )
  { 
    return const_cast<G4ParticleDefinition*>(ion);
  }
  else
  {
    return nullptr;
  }
}

// --------------------------------------------------------------------
//
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4int lvl)
{
  if ( (A<1) || (Z<=0) || (lvl<0) || (A>999) )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::FindIon(): illegal atomic number/mass"
             << " or excitation level:" << G4endl 
             << " Z =" << Z << "  A = " << A <<  "  IsoLvl = " << lvl << G4endl;
    }
#endif
    G4Exception("G4IonTable::FindIon()", "PART107",
                JustWarning, "illegal atomic number/mass");
    return nullptr;
  }
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion = nullptr;
  G4bool isFound = false;

  // check if light ion
  ion = GetLightIon(Z,A);
  if (ion != nullptr && lvl==0)
  { 
    // light ion 
    isFound = true;
  }
  else
  {    
    // -- loop over all particles in Ion table
    G4int encoding=GetNucleusEncoding(Z, A);
    for(auto i = fIonList->find(encoding); i != fIonList->cend(); ++i)
    {
      ion = i->second;
      if ( (ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
      // excitation level
      if ( ((const G4Ions*)(ion))->GetIsomerLevel() == lvl)
      {
        isFound = true;
        break;
      }
    } 
  }

  if ( isFound )
  { 
    if(lvl==9)
    {
      G4Exception("G4IonTable::FindIon()","PART5107", JustWarning,
                  "Isomer level 9 may be ambiguous.");
    }
    return const_cast<G4ParticleDefinition*>(ion);
  }
  else
  {
    return nullptr;
  }
}

// --------------------------------------------------------------------
//
G4ParticleDefinition*
G4IonTable::FindIon(G4int Z, G4int A, G4int LL, G4int lvl)
{
  if (LL==0) return FindIon(Z,A,lvl);
  
  if (A < 2 || Z < 0 || Z > A-LL || LL>A || A>999 )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::FindIon(): illegal atomic number/mass"
             << " or excitation level:" << G4endl
             << " Z =" << Z << "  A = " << A << " L = " << LL 
             <<"  IsomerLvl = " << lvl << G4endl;
    }    
#endif
    G4Exception( "G4IonTable::FindIon()", "PART107",
                 JustWarning, "illegal atomic number/mass");
    return nullptr;
  }

  // Search ions with A, Z ,E, lvl
  const G4ParticleDefinition* ion = nullptr;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A, LL);
  for(auto i = fIonList->find(encoding); i != fIonList->cend() ; ++i)
  {
    ion = i->second;
    if ( (ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if ( ion->GetQuarkContent(3) != LL) break;
    // excitation level
    if ( ((const G4Ions*)(ion))->GetIsomerLevel() == lvl)
    {
      isFound = true;
      break;
    }
  }

  if ( isFound )
  { 
    if(lvl==9)
    {
      G4Exception("G4IonTable::FindIon()", "PART5107", JustWarning,
                  "Isomer level 9 may be ambiguous.");
    }
    return const_cast<G4ParticleDefinition*>(ion);
  }
  else
  {
    return nullptr;
  }
}

// --------------------------------------------------------------------
// GetNucleusEncoding
//
G4int G4IonTable::GetNucleusEncoding(G4int Z, G4int A, G4double E, G4int lvl)
{
  // PDG code for Ions
  // Nuclear codes are given as 10-digit numbers +-100ZZZAAAI.
  // For a nucleus consisting of np protons and nn neutrons
  // A = np + nn and Z = np.
  // I gives the isomer level, with I = 0 corresponding 
  // to the ground state and I >0 to excitations
  
  if ( Z==1 && A==1 && E==0.0 ) return 2212; // proton
  
  G4int encoding = 1000000000;
  encoding += Z * 10000;
  encoding += A *10;
  if (lvl>0&&lvl<10) encoding +=lvl;       //isomer level
  else if (E>0.0) encoding += 9;  //isomer level
  
  return encoding;
}

// --------------------------------------------------------------------
// GetNucleusEncoding
//
G4int G4IonTable::GetNucleusEncoding(G4int Z, G4int A, G4int LL,
                                     G4double E, G4int lvl)
{
  // Get PDG code for Hyper-Nucleus Ions 
  // Nuclear codes are given as 10-digit numbers +-10LZZZAAAI.
  // For a nucleus consisting of np protons and nn neutrons
  // A = np + nn +nlambda and Z = np.
  // LL = nlambda
  // I gives the isomer level, with I = 0 corresponding 
  // to the ground state and I >0 to excitations

  G4int encoding = GetNucleusEncoding(Z, A, E, lvl);
  if (LL==0) return encoding; 
  encoding += LL*  10000000;
  if ( Z==1 && A==1 && E==0.0 ) encoding = 3122; // Lambda 

  return encoding;
}

// --------------------------------------------------------------------
// GetNucleusByEncoding
//
G4bool G4IonTable::GetNucleusByEncoding(G4int encoding,
                                        G4int& Z,    G4int& A, 
                                        G4double& E, G4int& lvl)
{
  if (encoding <= 0) return false; // anti particle   

  if (encoding == 2212)   // proton
  {
    Z = 1; A = 1;
    E = 0.0; lvl =0; 
    return true;
  } 

  encoding -= 1000000000;
  Z = encoding/10000;
  encoding -= 10000*Z;
  A = encoding/10;
  lvl = encoding % 10;
  return true;
}

// --------------------------------------------------------------------
// GetNucleusByEncoding
//
G4bool G4IonTable::GetNucleusByEncoding(G4int encoding,
                                        G4int& Z,    G4int& A, 
                                        G4int& LL,   
                                        G4double& E, G4int& lvl)
{
  if (encoding <= 0) return false; // anti particle   

  if (encoding == 3122)   // Lambda
  {
    Z = 1; A = 1; LL = 1;
    E = 0.0; lvl =0; 
    return true;
  } 

  if (encoding % 10 != 0)
  {
    // !!!not supported for excitation states !!!   
    return false;
  }
  if (encoding < 1000000000)
  {
    // anti particle   
    return false;
  }

  encoding -= 1000000000;
  LL = encoding/10000000;
  encoding -= 10000000*LL;
  Z = encoding/10000;
  encoding -= 10000*Z;
  A = encoding/10;
  lvl = encoding % 10;
  return true;
}

// --------------------------------------------------------------------
// GetIonName
//
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4double E,
                G4Ions::G4FloatLevelBase flb) const 
{
  static G4ThreadLocal G4String* pname = nullptr;
  if ( pname == nullptr )
  {
    pname = new G4String("");
    G4AutoDelete::Register(pname);
  }
  G4String& name = *pname;

  static G4ThreadLocal std::ostringstream* os = nullptr;
  if ( os == nullptr )
  {
    os = new std::ostringstream();
    G4AutoDelete::Register(os); 
    os->setf(std::ios::fixed);
    os->precision(3);
  }

  name = GetIonName(Z, A);

  // Excited energy
  if ( E>0  || flb!=G4Ions::G4FloatLevelBase::no_Float)
  {
    os->str("");
    std::ostringstream& oo = *os;

    // Excited nucleus
    oo<<'['<<E/keV;
    if (flb!=G4Ions::G4FloatLevelBase::no_Float)
    {
      oo<<G4Ions::FloatLevelBaseChar(flb);
    }
    oo<< ']';
    name += os->str();
  }

  return name;
}

// --------------------------------------------------------------------
// GetIonName
//
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4int LL, G4double E,
                G4Ions::G4FloatLevelBase flb) const 
{
  if (LL==0) return GetIonName(Z, A, E, flb); 
  static G4ThreadLocal G4String* pname = nullptr;
  if (pname == nullptr)
  {
    pname = new G4String("");
    G4AutoDelete::Register(pname);
  }
  G4String& name = *pname;
  name = "";
  for (G4int i=0; i<LL; ++i)
  {
    name +="L";
  }
  name += GetIonName(Z, A, E, flb);
  return name;
}

// --------------------------------------------------------------------
// GetIonName
//
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4int lvl) const 
{
  static G4ThreadLocal G4String* pname = nullptr;
  if ( pname == nullptr )
  {
    pname = new G4String("");
    G4AutoDelete::Register(pname);
  }
  G4String& name = *pname;

  static G4ThreadLocal std::ostringstream* os = nullptr;
  if ( os == nullptr )
  {
    os = new std::ostringstream();
    G4AutoDelete::Register(os);
    os->setf(std::ios::fixed);
  }

  if ( (0< Z) && (Z <=numberOfElements) )
  {
    name = elementName[Z-1];
  }
  else if (Z > numberOfElements)
  {
    os->str("");
    os->operator<<(Z);
    name = "E" + os->str() + "-";
  }
  else
  {
    name = "?";
    return name;
  }
  // Atomic Mass
  os->str("");
  os->operator<<(A);

  if ( lvl>0 )
  {
    std::ostringstream& oo = *os;
    // Isomer level for Excited nucelus 
    oo<<'['<<lvl << ']';
  }
  name += os->str();

  return name;
}

// --------------------------------------------------------------------
// GetIonName
//
const G4String&
G4IonTable::GetIonName(G4int Z, G4int A, G4int LL, G4int lvl) const 
{
  if (LL==0) return GetIonName(Z, A, lvl); 
  static G4ThreadLocal G4String* pname = nullptr;
  if ( pname == nullptr )
  {
    pname = new G4String("");
    G4AutoDelete::Register(pname);
  }
  G4String &name = *pname;
  for (G4int i=0; i<LL; ++i)
  {
    name +="L";
  }
  name += GetIonName(Z, A, lvl);
  return name;
}

// --------------------------------------------------------------------
// IsIon
//
G4bool G4IonTable::IsIon(const G4ParticleDefinition* particle)
{
  // Return true if the particle is ion
  static const G4String nucleus("nucleus");
  static const G4String proton("proton");

  // Neutron is not ion
  if ( (particle->GetAtomicMass()>0)
    && (particle->GetAtomicNumber()>0) )
  {
    if (particle->GetBaryonNumber()>0)  return true;
    else return false;
  }
   
  // Particles derived from G4Ions
  if (particle->GetParticleType() == nucleus) return true;

  // Proton (Hydrogen nucleus)
  if (particle->GetParticleName() == proton) return true;

  return false;
}

// --------------------------------------------------------------------
// IsAntiIon
//
G4bool G4IonTable::IsAntiIon(const G4ParticleDefinition* particle)
{
  // Return true if the particle is ion
  static const G4String anti_nucleus("anti_nucleus");
  static const G4String anti_proton("anti_proton");

  // Anti_neutron is not ion
  if ( (particle->GetAtomicMass()>0)
    && (particle->GetAtomicNumber()>0) )
  {
    if (particle->GetBaryonNumber()<0)  return true;
    else return false;
  }

  // Particles derived from G4Ions
  if (particle->GetParticleType() == anti_nucleus) return true;

  // Anti_proton (Anti_Hydrogen nucleus)
  if (particle->GetParticleName() == anti_proton) return true;

  return false;
}

// --------------------------------------------------------------------
// IsLightIon
//
G4bool G4IonTable::IsLightIon(const G4ParticleDefinition* particle) const
{
  static const std::string names[]
  = { "proton", "alpha", "deuteron", "triton", "He3"};

  // Return true if the particle is pre-defined ion
  return std::find(names, names+5, (particle->GetParticleName()).c_str())!=names+5;
} 

// --------------------------------------------------------------------
// IsLightAntiIon
//
G4bool G4IonTable::IsLightAntiIon(const G4ParticleDefinition* particle) const
{
  static const std::string names[]
  = { "anti_proton", "anti_alpha", "anti_deuteron", "anti_triton", "anti_He3"};

  // Return true if the particle is pre-defined ion
  return std::find(names, names+5, (particle->GetParticleName()).c_str())!=names+5;
} 

// --------------------------------------------------------------------
// GetLightIon
//
G4ParticleDefinition* G4IonTable::GetLightIon(G4int Z, G4int A) const
{
  // Returns pointer to pre-defined ions
  const G4ParticleDefinition* ion = nullptr;
  if ( (Z<=2) )
  {
#ifndef G4MULTITHREADED
    // In sequential use lazy-initialization
    lightions::Init();
#endif
    if ( (Z==1)&&(A==1) ) {
      ion = lightions::p_proton;
    } else if ( (Z==1)&&(A==2) ) {
        ion = lightions::p_deuteron;
    } else if ( (Z==1)&&(A==3) ) {
      ion = lightions::p_triton;
    } else if ( (Z==2)&&(A==4) ) {
      ion = lightions::p_alpha;
    } else if ( (Z==2)&&(A==3) ) {
      ion = lightions::p_He3;
    }
  }
  return const_cast<G4ParticleDefinition*>(ion);
}

// --------------------------------------------------------------------
// GetLightAntiIon
//
G4ParticleDefinition* G4IonTable::GetLightAntiIon(G4int Z, G4int A) const
{
  // Returns pointer to pre-defined ions 
  const G4ParticleDefinition* ion = nullptr;
  if ( (Z<=2) )
  {
#ifndef G4MULTITHREADED
    // In sequential use lazy-initialization
    antilightions::Init();
#endif
    if ( (Z==1)&&(A==1) ) {
      ion = antilightions::p_proton;
    } else if ( (Z==1)&&(A==2) ) {
      ion = antilightions::p_deuteron;
    } else if ( (Z==1)&&(A==3) ) {
      ion = antilightions::p_triton;
    } else if ( (Z==2)&&(A==4) ) {
      ion = antilightions::p_alpha;
    } else if ( (Z==2)&&(A==3) ) {
      ion = antilightions::p_He3;
    }
  }
  return const_cast<G4ParticleDefinition*>(ion);
}

// --------------------------------------------------------------------
// GetNucleusMass
//
G4double
G4IonTable::GetNucleusMass(G4int Z, G4int A, G4int LL, G4int lvl) const
{
  if ( (A<1)  || (Z<0) || (LL<0) || (lvl<0) || (lvl>9) )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4IonTable::GetNucleusMass() : illegal atomic number/mass:"
             << G4endl
             << " Z =" << Z << "  A = " << A  
             << " L = " << LL << " lvl = " << lvl << G4endl;
    }
#endif
    G4Exception("G4IonTable::GetNucleusMass()","PART107",
                EventMustBeAborted, "illegal atomic number/mass");
    return -1.0;
  }
  
  G4double mass;
  if (LL == 0)
  {
    // calculate nucleus mass
    const G4ParticleDefinition* ion=GetLightIon(Z, A);
    
    if (ion != nullptr)
    {
      mass = ion->GetPDGMass();
    }
    else
    {
      // Use G4NucleiProperties::GetNuclearMass
      mass = G4NucleiProperties::GetNuclearMass(A, Z);
    }
    
    // Isomer
    if ( lvl>0 )
    {
      // -- loop over all particles in Ion table
      G4int encoding=GetNucleusEncoding(Z, A);
      G4bool isFound = false;
      for(auto i = fIonList->find(encoding);i != fIonList->cend() ; ++i)
      {
        ion = i->second;
        if ( ( ion->GetAtomicNumber()!=Z) || (ion->GetAtomicMass()!=A) ) break;
        // Excitation level
        if ( ((const G4Ions*)(ion))->GetIsomerLevel() == lvl)
        {
          isFound = true;
          break;
        }
      }
      if (isFound)
      {
        // Return existing isomer mass
        mass = ion->GetPDGMass();
      }
      else
      {
        // Find isomer from IsotopeTable
        const G4IsotopeProperty*  fProperty = FindIsotope(Z, A, lvl);
        if (fProperty != nullptr ) mass += fProperty->GetEnergy();
      }
    }
  }
  else
  {
    mass = G4HyperNucleiProperties::GetNuclearMass(A, Z, LL);
  }
  return mass;
}

// --------------------------------------------------------------------
// GetIsomerMass
//
G4double G4IonTable::GetIsomerMass(G4int Z, G4int A, G4int  lvl) const
{
  return GetNucleusMass(Z,A,0,lvl);
}

// --------------------------------------------------------------------
// GetIonMass
//
G4double G4IonTable::GetIonMass(G4int Z, G4int A, G4int LL, G4int lvl) const
{
  return GetNucleusMass(Z,A,LL,lvl);
}

// --------------------------------------------------------------------
// -- Methods for handling container  ---
// --------------------------------------------------------------------
//
void G4IonTable::clear()
{
  if (G4ParticleTable::GetParticleTable()->GetReadiness())
  {
    G4Exception("G4IonTable::clear()",
                "PART116", JustWarning,
                "No effects because readyToUse is true.");    
    return;
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel()>2)
  {
    G4cout << "G4IonTable::Clear() : number of Ion registered =  "; 
    G4cout << fIonList->size() <<  G4endl;
  }
#endif
  fIonList->clear();
}

// --------------------------------------------------------------------
//
void G4IonTable::Insert(const G4ParticleDefinition* particle)
{
  if (!IsIon(particle)) return;
  if (Contains(particle)) return;
 
  G4int Z = particle->GetAtomicNumber();
  G4int A = particle->GetAtomicMass();  
  G4int LL = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, LL); // encoding of the groud state 

  // Register the ion with its encoding of the ground state  
  fIonListShadow->insert( std::pair<const G4int,
                          const G4ParticleDefinition*>(encoding, particle) );
}

// --------------------------------------------------------------------
//
void G4IonTable::InsertWorker(const G4ParticleDefinition* particle)
{
  if(!particle) return;

  G4int Z = particle->GetAtomicNumber();
  G4int A = particle->GetAtomicMass();  
  G4int LL = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, LL);
  G4bool found = false;
  if (encoding !=0 )
  {
    for(auto i = fIonList->find(encoding); i != fIonList->cend(); ++i)
    {
      if (particle == i->second)
      {
        found  = true;
        break;
      }
    }
  }
  if(found) return;
 
  // Register the ion with its encoding of the gronud state  
  fIonList->insert( std::pair<const G4int,
                    const G4ParticleDefinition*>(encoding, particle) );
}

// --------------------------------------------------------------------
//
void G4IonTable::Remove(const G4ParticleDefinition* particle)
{
  if(particle == nullptr) return;
#ifdef G4MULTITHREADED
  if(G4Threading::IsWorkerThread())
  {
    G4ExceptionDescription ed;
    ed << "Request of removing " << particle->GetParticleName()
       << " is ignored as it is invoked from a worker thread.";
    G4Exception("G4IonTable::Remove()", "PART10117", JustWarning, ed);
    return;
  }
#endif
  if (G4ParticleTable::GetParticleTable()->GetReadiness())
  {
    G4StateManager* pStateManager = G4StateManager::GetStateManager();
    G4ApplicationState currentState = pStateManager->GetCurrentState();
    if (currentState != G4State_PreInit)
    {
      G4String msg = "Request of removing ";
      msg += particle->GetParticleName();  
      msg += " has No effects other than Pre_Init";
      G4Exception("G4IonTable::Remove()",
                  "PART117", JustWarning, msg);
      return;
    }
    else
    {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
      {
        G4cout << particle->GetParticleName()
               << " will be removed from the IonTable " << G4endl;
      }
#endif
    }
  }

  if (IsIon(particle))
  {
    G4int Z = particle->GetAtomicNumber();
    G4int A = particle->GetAtomicMass();  
    G4int LL = particle->GetQuarkContent(3);  // strangeness 
    G4int encoding=GetNucleusEncoding(Z, A, LL);
    if (encoding !=0 )
    {
      for(auto i = fIonListShadow->find(encoding);
               i != fIonListShadow->cend() ; ++i)
      {
        if (particle == i->second)
        {
          fIonListShadow->erase(i);
          break;
        }
      }
    }
  }
  else
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1)
    {
      G4cout << "G4IonTable::Remove :" << particle->GetParticleName() 
             << " is not ions" << G4endl; 
    }
#endif
  }
}

// --------------------------------------------------------------------
// -- Dump Information 
//
void G4IonTable::DumpTable(const G4String& particle_name) const
{
  const G4ParticleDefinition* ion;
  for (auto idx = fIonList->cbegin(); idx!= fIonList->cend(); ++idx)
  {
    ion = idx->second;
    if (( particle_name == "ALL" ) || (particle_name == "all"))
    {
      ion->DumpTable();
    }
    else if ( particle_name == ion->GetParticleName() )
    {
      ion->DumpTable();
    }
  }
}

// --------------------------------------------------------------------
//
const G4String G4IonTable::elementName[] =
{
  "H",                                                                                                "He", 
  "Li", "Be",                                                             "B",  "C",  "N",  "O", "F", "Ne", 
  "Na", "Mg",                                                             "Al", "Si", "P", "S", "Cl", "Ar", 
  "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", 
  "Rb", "Sr", "Y", "Zr", "Nb", "Mo","Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", 
  "Cs", "Ba", 
              "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", 
                   "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", 
  "Fr", "Ra", 
              "Ac", "Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
              "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};

// --------------------------------------------------------------------
// GetVerboseLevel
//
G4int G4IonTable::GetVerboseLevel() const
{
  return G4ParticleTable::GetParticleTable()->GetVerboseLevel();
}

// --------------------------------------------------------------------
// AddProcessManager
//
void  G4IonTable::AddProcessManager(G4ParticleDefinition* ion)
{
  if(ion->IsGeneralIon())
  {
    // Check whether GenericIon has processes
    G4ParticleDefinition* genericIon = 
      G4ParticleTable::GetParticleTable()->GetGenericIon();

    G4ProcessManager* pman = nullptr;
    if (genericIon != nullptr) pman = genericIon->GetProcessManager();
    if ((genericIon == nullptr) || (genericIon->GetParticleDefinitionID() < 0)
                                || (pman==nullptr))
    {
      G4String msg = "G4IonTable::AddProcessManager(): cannot create ion of ";
      msg += ion->GetParticleName();
      msg += "\n because GenericIon is not available!!";
      G4Exception("G4IonTable::AddProcessManager()", "PART105",
                  FatalException, msg);
      return;
    }
  
    ion->SetParticleDefinitionID(genericIon->GetParticleDefinitionID());
  }
  else
  {
    // Is this a MuonicAtom ?
    G4MuonicAtom* muatom = dynamic_cast<G4MuonicAtom*> (ion);

    if ( muatom != nullptr )
    {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1)
      {
        G4cout << "G4IonTable::AddProcessManager(): "
               << "MuonicAtom dynamic_cast succeeded for " 
               << ion->GetParticleName() << G4endl;
      }
#endif
      // Check whether GenericMuonicAtom has processes
      G4ParticleDefinition* genericMA = 
        G4ParticleTable::GetParticleTable()->GetGenericMuonicAtom();

      G4ProcessManager* pman = nullptr;
      if (genericMA != nullptr) pman = genericMA->GetProcessManager();
      if ((genericMA == nullptr) || (genericMA->GetParticleDefinitionID() < 0)
                                 || (pman==nullptr))
      {
        G4String msg =
               "G4IonTable::AddProcessManager(): cannot create MuonicAtom ";
        msg += ion->GetParticleName();
        msg += "\n because GenericMuonicAtom is not available!!";
        G4Exception("G4IonTable::AddProcessManager()",
                    "PART106", FatalException, msg);
        return;
      }
  
      ion->SetParticleDefinitionID(genericMA->GetParticleDefinitionID());
    }
    else
    {
      G4String msg =
               "G4IonTable::AddProcessManager(): cannot create ";
      msg += ion->GetParticleName();
      msg += "\n because of unsupported particle type !!";
      G4Exception("G4IonTable::AddProcessManager()", "PART107",
                  FatalException, msg);
      return;
    }
  }
  return;
}

// --------------------------------------------------------------------
// RegisterIsotopeTable
//
void  G4IonTable::RegisterIsotopeTable(G4VIsotopeTable* table)
{
  //check duplication
  G4String name = table->GetName();
  for (std::size_t i=0; i<fIsotopeTableList->size(); ++i)
  {
    G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[i];
    if (name == fIsotopeTable->GetName()) return;
  }
  // register 
  fIsotopeTableList->push_back(table);
}

// --------------------------------------------------------------------
// GetIsotopeTable
//
G4VIsotopeTable* G4IonTable::GetIsotopeTable(std::size_t index) const
{
   G4VIsotopeTable* fIsotopeTable = nullptr;
   if ( index < fIsotopeTableList->size() )
   {
     fIsotopeTable = (*fIsotopeTableList)[index];
   }
   return fIsotopeTable;
}

// --------------------------------------------------------------------
// FindIsotope
//
G4IsotopeProperty* G4IonTable::FindIsotope(G4int Z, G4int A, G4double E,
                                           G4Ions::G4FloatLevelBase flb) const
{
  if (fIsotopeTableList == nullptr) return nullptr;
  if (fIsotopeTableList->size() == 0) return nullptr;
  
  G4IsotopeProperty* property = nullptr;

  for (std::size_t i=0; i<fIsotopeTableList->size(); ++i)
  {
    G4VIsotopeTable* fIsotopeTable
      = (*fIsotopeTableList)[fIsotopeTableList->size()-i-1];
    property = fIsotopeTable->GetIsotope(Z,A,E,flb);
    if(property) break;
  }
  
  return property;
}

// --------------------------------------------------------------------
// FindIsotope
//
G4IsotopeProperty* G4IonTable::FindIsotope(G4int Z, G4int A, G4int lvl) const
{
  if (fIsotopeTableList == nullptr) return nullptr;
  if (fIsotopeTableList->size()==0) return nullptr;
  
  G4IsotopeProperty* property = nullptr;

  // iterate 
  for (std::size_t i=0; i<fIsotopeTableList->size(); ++i)
  {
    G4VIsotopeTable* fIsotopeTable
      = (*fIsotopeTableList)[fIsotopeTableList->size()-i-1];
    property = fIsotopeTable->GetIsotope(Z,A,lvl);
    if(property) break;
  }
  
  return property;
}

// --------------------------------------------------------------------
// CreateAllIon
//
void G4IonTable::CreateAllIon()
{
  PreloadNuclide();
}

// --------------------------------------------------------------------
// CreateAllIsomer
//
void G4IonTable::CreateAllIsomer()
{
  PreloadNuclide();
}

// --------------------------------------------------------------------
// PrepareNuclideTable
//
void G4IonTable::PrepareNuclideTable()
{
  if (pNuclideTable == nullptr)
    pNuclideTable = G4NuclideTable::GetNuclideTable();
}

// --------------------------------------------------------------------
// PreloadNuclide
//
void G4IonTable::PreloadNuclide()
{
  if ( isIsomerCreated || !G4Threading::IsMultithreadedApplication() ) return;

  pNuclideTable->GenerateNuclide();

  for ( std::size_t i=0 ; i!=pNuclideTable->entries(); ++i )
  {
    const G4IsotopeProperty* fProperty = pNuclideTable->GetIsotopeByIndex( i );
    G4int Z = fProperty->GetAtomicNumber();
    G4int A = fProperty->GetAtomicMass();
    G4double Eex = fProperty->GetEnergy();
    GetIon(Z,A,Eex);
  }

  isIsomerCreated = true;
}

// --------------------------------------------------------------------
// GetParticle
//
G4ParticleDefinition* G4IonTable::GetParticle(G4int index) const
{
  if ( (index >=0) && (index < Entries()) )
  {
    auto idx = fIonList->cbegin();
    G4int counter = 0;
    while( idx != fIonList->cend() ) // Loop checking, 09.08.2015, K.Kurashige
    {
      if ( counter == index )
      {
        return const_cast<G4ParticleDefinition*>(idx->second);
      }
      ++counter;
      ++idx;
    }
  } 
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
  {
    G4cout << " G4IonTable::GetParticle"
           << " invalid index (=" << index << ")" 
           << " entries = " << Entries() << G4endl;
  }
#endif
  return nullptr;
}

// --------------------------------------------------------------------
// Contains
//
G4bool G4IonTable::Contains(const G4ParticleDefinition* particle) const
{
  if (!IsIon(particle)) return false;

  G4int Z = particle->GetAtomicNumber();
  G4int A = particle->GetAtomicMass();  
  G4int LL = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, LL);
  G4bool found = false;
  if (encoding != 0 )
  {
    for(auto i = fIonListShadow->find(encoding);
             i != fIonListShadow->cend(); ++i)
    {
      if (particle == i->second )
      {
        found  = true;
        break;
      }
    }
  }
  return found;
}

// --------------------------------------------------------------------
// Entries
//
G4int G4IonTable::Entries() const
{
  return (G4int)fIonList->size();
}

// --------------------------------------------------------------------
// size
//
G4int G4IonTable::size() const
{
  return (G4int)fIonList->size();
}

// --------------------------------------------------------------------
// FindIonInMaster
//
G4ParticleDefinition*
G4IonTable::FindIonInMaster(G4int Z, G4int A, G4double E, 
                            G4Ions::G4FloatLevelBase flb, G4int /*J*/)
{
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion = nullptr;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A);
  for(auto i = fIonListShadow->find(encoding); i != fIonListShadow->cend(); ++i)
  {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    // excitation level
    G4double anExcitaionEnergy = ((const G4Ions*)(ion))->GetExcitationEnergy();
    if (std::fabs(E - anExcitaionEnergy) < pNuclideTable->GetLevelTolerance() )
    {
      if(((const G4Ions*)(ion))->GetFloatLevelBase()==flb)
      {
        isFound = true;
        break;
      }
    }
  }

  if ( isFound )
  { 
    return const_cast<G4ParticleDefinition*>(ion);
  }
  else
  {
    return nullptr;
  }
}

// --------------------------------------------------------------------
// FindIonInMaster
//
G4ParticleDefinition*
G4IonTable::FindIonInMaster(G4int Z, G4int A, G4int LL, G4double E,
                            G4Ions::G4FloatLevelBase flb, G4int J)
{
  if (LL==0) return FindIon(Z,A,E,flb,J);
  
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion = nullptr;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding = GetNucleusEncoding(Z, A, LL, 0.0, 0);
  for(auto i = fIonListShadow->find(encoding); i != fIonListShadow->cend(); ++i)
  {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if(  ion->GetQuarkContent(3) != LL) break;
    // Excitation level
    G4double anExcitaionEnergy = ((const G4Ions*)(ion))->GetExcitationEnergy();
    if (std::fabs(E - anExcitaionEnergy) < pNuclideTable->GetLevelTolerance() )
    {
      if(((const G4Ions*)(ion))->GetFloatLevelBase()==flb)
      {
        isFound = true;
        break;
      }
    }
  }

  if ( isFound )
  { 
    return const_cast<G4ParticleDefinition*>(ion);
  }
  else
  {
    return nullptr;
  }
}

// --------------------------------------------------------------------
// FindIonInMaster
//
G4ParticleDefinition* G4IonTable::FindIonInMaster(G4int Z, G4int A, G4int lvl)
{
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion = nullptr;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A);
  for(auto i = fIonListShadow->find(encoding); i != fIonListShadow->cend(); ++i)
  {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    // Excitation level
    if ( ((const G4Ions*)(ion))->GetIsomerLevel() == lvl)
    {
      isFound = true;
      break;
    } 
  }

  if ( isFound )
  { 
    return const_cast<G4ParticleDefinition*>(ion);
  }
  else
  {
    return nullptr;
  }
}

// --------------------------------------------------------------------
// FindIonInMaster
//
G4ParticleDefinition*
G4IonTable::FindIonInMaster(G4int Z, G4int A, G4int LL, G4int lvl)
{
  if (LL==0) return FindIon(Z,A,lvl);
  
  // Search ions with A, Z ,E, lvl
  const G4ParticleDefinition* ion = nullptr;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A, LL);
  for(auto i = fIonListShadow->find(encoding); i != fIonListShadow->cend(); ++i)
  {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if ( ion->GetQuarkContent(3) != LL) break;
    // excitation level
    if ( ((const G4Ions*)(ion))->GetIsomerLevel() == lvl)
    {
      isFound = true;
      break;
    }
  }

  if ( isFound )
  { 
    return const_cast<G4ParticleDefinition*>(ion);
  }
  else
  {
    return nullptr;
  }
}

// --------------------------------------------------------------------
// GetLifeTime
//
G4double G4IonTable::GetLifeTime(const G4ParticleDefinition* particle) const
{
  if((particle->IsGeneralIon()) && (pNuclideTable == nullptr))
  {
   G4Exception("G4IonTable::GetLifeTime()", "ParticleIon1001", FatalException,
               "Method is invoked before G4IonTable is initialized.");
   return 0.;
  } 
  return particle->GetPDGLifeTime();
}

// --------------------------------------------------------------------
// GetLifeTime
//
G4double
G4IonTable::GetLifeTime(G4int Z, G4int A, G4double E, char flbChar) const
{ 
  return GetLifeTime(Z,A,E,G4Ions::FloatLevelBase(flbChar)); 
}

// --------------------------------------------------------------------
// GetLifeTime
//
G4double G4IonTable::GetLifeTime(G4int Z, G4int A, G4double E,
             G4Ions::G4FloatLevelBase flb) const
{
  G4double life = -1001.0;
  const G4IsotopeProperty* fProperty = FindIsotope(Z, A, E, flb);
  if( fProperty != nullptr ) life = fProperty->GetLifeTime();
  return life;
}

// --------------------------------------------------------------------
// GetMuonicAtom
//
G4ParticleDefinition* G4IonTable::GetMuonicAtom(G4Ions const* base)
{
  if (base==nullptr || !IsIon(base))
  {
    G4Exception("G4IonTable::GetMuonicAtom()", "PART987654321",
                FatalException, "Constructor argument is not a G4Ions");
    return nullptr;
  }

  // We're assuming here that we get a base that is actually
  // constructed and unexcited ... strip excitations, Lambdas, and
  // isomers from the encoding

  auto const Z = base->GetAtomicNumber();
  auto const A = base->GetAtomicMass();
  auto const baseenc = GetNucleusEncoding(Z,A);
  auto const encoding = baseenc+1000000000;

  // We have to do all the MT manipulations manually, because the
  // convenience functions assume a G4Ions with canonical PDG codes;
  // they recalculate the encoding from particle properties rather
  // than using the carried member function values.  Thus, they will
  // do operations on the base ion, rather than the passed in
  // G4MuonicAtom

  auto i = fIonList->find(encoding);
  if(i!=fIonList->cend())
  {
    return const_cast<G4ParticleDefinition*>(i->second);
  }
  // not in threadlocal list; check global list ... 
#ifdef G4MULTITHREADED
  if(G4Threading::IsWorkerThread())
  {
    G4MUTEXLOCK(&G4IonTable::ionTableMutex);
    i = fIonListShadow->find(encoding);
    auto end = fIonListShadow->cend();
    G4MUTEXUNLOCK(&G4IonTable::ionTableMutex);
    if(i!=end)
    {
      // we found it, stuff it into the threadlocal list
      fIonList->insert(*i);
      // and then return it ... 
      return const_cast<G4ParticleDefinition*>(i->second);
    }
  }
#endif

  // not found in either list; create and potentially insert 
  auto const name = "Mu"+GetIonName(Z,A);

  G4MuonicAtom* muatom = 
    G4MuonicAtomHelper::ConstructMuonicAtom(name, encoding, base);

  // Not sure this is doing the right thing...
  AddProcessManager(muatom);

  // Now, we have to push the muatom into the appropriate IonTables
  // first, recheck global list, in case another thread came along
  // before us and created this same muatom

#ifdef G4MULTITHREADED
  if(G4Threading::IsWorkerThread())
  {
    G4MUTEXLOCK(&G4IonTable::ionTableMutex);
    // first, we need to make sure it hasn't been inserted by some
    // other thread
    auto j = fIonListShadow->find(encoding);
    if( j!= fIonListShadow->cend() )
    {
      // oops ... someone else built a copy when we weren't looking;
      // cleanup our instantiation, and take a handle to the one in
      // the global list
      delete muatom;
      muatom = const_cast<G4MuonicAtom*>
               (static_cast<G4MuonicAtom const*>(j->second));
    }
    else
    {
      // otherwise, push onto the global list first
      fIonListShadow->insert(std::make_pair(encoding, muatom));
    }
    G4MUTEXUNLOCK(&G4IonTable::ionTableMutex);
  }
#endif  
  // in either case, push onto the the threadlocal list
  fIonList->insert(std::make_pair(encoding,muatom));

  return muatom; 
}

// --------------------------------------------------------------------
// GetMuonicAtom
//
G4ParticleDefinition* G4IonTable::GetMuonicAtom(G4int Z, G4int A)
{
  // Need the cast because we need a G4Ions* to pass into the
  // function, but GetIon returns a G4ParticleDefinition* 
  auto base = static_cast<G4Ions const*>(GetIon(Z,A, 0.0));
  return GetMuonicAtom(base);
}
