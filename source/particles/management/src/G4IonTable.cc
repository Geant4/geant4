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
//
// $Id: G4IonTable.cc 102916 2017-03-02 12:58:58Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	History: first implementation, based on object model of
//	27 June 1998  H.Kurashige
// ---------------------------------------------------------------
//      modified GetIon                 02 Aug., 98 H.Kurashige
//      added Remove()                  06 Nov.,98 H.Kurashige
//      use G4NucleiPropoerties to get nuceli Mass 17  Nov.,98 H.Kurashige
//      use G4GenericIon for process List
//      modify fomula of Ion mass       09 Dec., 98 H.Kurashige 
//          -----
//      Modified GetIon methods         17 Aug. 99 H.Kurashige
//      New design using G4VIsotopeTable          5 Oct. 99 H.Kurashige
//      Modified Element Name for Z>103  06 Apr. 01 H.Kurashige
//      Remove test of cuts in SetCuts   16 Jan  03 V.Ivanchenko
//      Added initial support for Muonic Atoms   1 Jul 16  K.Lynch

#include <iostream>               
#include <iomanip>               
#include <sstream>

#include "G4ios.hh"
#include "G4Threading.hh"

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
G4ThreadLocal G4IonTable::G4IonList* G4IonTable::fIonList = 0;
G4ThreadLocal std::vector<G4VIsotopeTable*> *G4IonTable::fIsotopeTableList = 0;
G4IonTable::G4IonList* G4IonTable::fIonListShadow = 0;
std::vector<G4VIsotopeTable*> *G4IonTable::fIsotopeTableListShadow = 0;

namespace lightions {
  static const G4ParticleDefinition* p_proton=0;
  static const G4ParticleDefinition* p_deuteron=0;
  static const G4ParticleDefinition* p_triton=0;
  static const G4ParticleDefinition* p_alpha=0;
  static const G4ParticleDefinition* p_He3=0;
  void Init() {
    if ( p_proton ) return;
    p_proton   = G4ParticleTable::GetParticleTable()-> FindParticle("proton"); // proton
    p_deuteron = G4ParticleTable::GetParticleTable()-> FindParticle("deuteron"); // deuteron
    p_triton   = G4ParticleTable::GetParticleTable()-> FindParticle("triton"); // tritoon
    p_alpha    = G4ParticleTable::GetParticleTable()-> FindParticle("alpha"); // alpha
    p_He3      = G4ParticleTable::GetParticleTable()-> FindParticle("He3"); // He3
  }
}

namespace antilightions {
    static const G4ParticleDefinition* p_proton=0;
    static const G4ParticleDefinition* p_deuteron=0;
    static const G4ParticleDefinition* p_triton=0;
    static const G4ParticleDefinition* p_alpha=0;
    static const G4ParticleDefinition* p_He3=0;
    void Init() {
        if ( p_proton ) return;
        p_proton   = G4ParticleTable::GetParticleTable()-> FindParticle("anti_proton"); // proton
        p_deuteron = G4ParticleTable::GetParticleTable()-> FindParticle("anti_deuteron"); // deuteron
        p_triton   = G4ParticleTable::GetParticleTable()-> FindParticle("anti_triton"); // tritoon
        p_alpha    = G4ParticleTable::GetParticleTable()-> FindParticle("anti_alpha"); // alpha
        p_He3      = G4ParticleTable::GetParticleTable()-> FindParticle("anti_He3"); // He3
    }
}

#ifdef G4MULTITHREADED
G4Mutex G4IonTable::ionTableMutex = G4MUTEX_INITIALIZER;
#endif 

////////////////////
G4IonTable::G4IonTable()
  : pNuclideTable(0),
    isIsomerCreated(false),
    n_error(0)
{
  fIonList = new G4IonList();

  // Set up the shadow pointer used by worker threads.
  //
  if (fIonListShadow == 0)
  {
    fIonListShadow = fIonList;
  }

  fIsotopeTableList = new std::vector<G4VIsotopeTable*>;

  // Set up the shadow pointer used by worker threads.
  //
  if (fIsotopeTableListShadow == 0)
  {
    fIsotopeTableListShadow = fIsotopeTableList;
  }    

  PrepareNuclideTable();
  RegisterIsotopeTable(pNuclideTable);
}

// This method is used by each worker thread to copy the content
// from the master thread.
//
void G4IonTable::SlaveG4IonTable()
{
G4Exception("G4IonTable::SlaveG4ParticleTable()","G4MT0000",FatalException,"Obsolete");
}

void G4IonTable::WorkerG4IonTable()
{
  if( fIonList == 0 )
  { fIonList = new G4IonList(); }
  else
  { fIonList->clear(); }

  G4IonListIterator it;
  for (it = fIonListShadow->begin() ; it != fIonListShadow->end(); it++ ) {
///////////////////////////    G4ParticleDefinition* ion = const_cast<G4ParticleDefinition*>(it->second);
///////////////////////////    if (ion->IsGeneralIon()) AddProcessManager(ion); 
    fIonList->insert(*it);
  }

  // Do not copy Isotoper Table to Worker thread
  if( fIsotopeTableList == 0 ) {
    fIsotopeTableList = new std::vector<G4VIsotopeTable*>; 
    for (size_t i = 0; i < fIsotopeTableListShadow->size(); i++){
      fIsotopeTableList->push_back((*fIsotopeTableListShadow)[i]); 
    }
  }

  /////////fIsotopeTableList = new std::vector<G4VIsotopeTable*>;
  /////////RegisterIsotopeTable(pNuclideTable);
}

void G4IonTable::InitializeLightIons()
{
  lightions::Init();
  antilightions::Init();
}


////////////////////
G4IonTable::~G4IonTable()
{
  // delete IsotopeTable if exists
  if (fIsotopeTableList != 0)
  {
    for (size_t i = 0; i< fIsotopeTableList->size(); ++i)
    {
      G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[i];
      //delete fIsotopeTable;
      if( fIsotopeTable != G4NuclideTable::GetNuclideTable() ) delete fIsotopeTable;
    }
    fIsotopeTableList->clear();
    delete fIsotopeTableList;
  }
  fIsotopeTableList =0;


  if (fIonList ==0) return;
  // remove all contents in the Ion List 
  //  No need to delete here because all particles are dynamic objects
  fIonList->clear();

  delete fIonList;
  fIonList =0;
}

////////////////////
void G4IonTable::DestroyWorkerG4IonTable()
{
  // delete IsotopeTable if exists
  if (fIsotopeTableList != 0)
  {
    for (size_t i = 0; i< fIsotopeTableList->size(); ++i)
    {
      G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[i];
      //delete fIsotopeTable;
      if( fIsotopeTable != G4NuclideTable::GetNuclideTable() ) delete fIsotopeTable;
    }
    fIsotopeTableList->clear();
    delete fIsotopeTableList;
  }
  fIsotopeTableList =0;


  if (fIonList ==0) return;
  // remove all contents in the Ion List 
  //  No need to delete here because all particles are dynamic objects
  fIonList->clear();

  delete fIonList;
  fIonList =0;
}


////////////////////
// -- CreateIon method ------
////////////////////
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, G4double E, 
                                            G4Ions::G4FloatLevelBase flb)
{
  G4ParticleDefinition* ion=0;

  // check whether GenericIon has processes
  G4ParticleDefinition* genericIon = 
    G4ParticleTable::GetParticleTable()->GetGenericIon();
  G4ProcessManager* pman=0;
  if (genericIon!=0) pman = genericIon->GetProcessManager();
  if ((genericIon ==0) || (genericIon->GetParticleDefinitionID() < 0) || (pman==0)){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "G4IonTable::CreateIon() : can not create ion of  " 
             << " Z =" << Z << "  A = " << A 
             << "  because GenericIon is not ready !!" <<   G4endl;
    }
#endif
    G4Exception( "G4IonTable::CreateIon()","PART105",
		 JustWarning, 
		 "Can not create ions because GenericIon is not ready");
    return 0;
  }
  
  G4double life = 0.0;
  G4DecayTable* decayTable =0;
  G4bool stable = true;
  G4double mu = 0.0;
  G4double Eex = 0.0;
  G4int    lvl =0;
  G4int    J=0;

  const G4IsotopeProperty*  fProperty = FindIsotope(Z, A, E, flb);
  if (fProperty !=0 ){
    Eex  = fProperty->GetEnergy();
    lvl  = fProperty->GetIsomerLevel();
    J    = fProperty->GetiSpin();
    life = fProperty->GetLifeTime();
    mu   = fProperty->GetMagneticMoment();    
    decayTable = fProperty->GetDecayTable();
    stable = (life <= 0.) || (decayTable ==0);
    lvl = fProperty->GetIsomerLevel();
    if (lvl <0) lvl=9;
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4ExceptionDescription ed;
      ed << "G4IonTable::CreateIon() : G4IsotopeProperty object was not found for"
         << " Z = " << Z << " A = " << A << " E = " << E/keV << " (keV)";
      if(flb!=G4Ions::G4FloatLevelBase::no_Float)
      { ed << " FloatingLevel +" << G4Ions::FloatLevelBaseChar(flb); }
      ed << ".\n"
         << " Physics quantities such as life are not set for this ion.";
      G4Exception( "G4IonTable::CreateIon()","PART70105", JustWarning, ed);
    }
#endif
    // excitation energy
    Eex = E;
    // lvl is assigned to 9 temporally    
    if (Eex>0.0) lvl=9;
  }

  //Eex = G4NuclideTable::Round(Eex); 
  if (Eex==0.0) lvl=0;
  // ion name
  G4String name =""; 
  /////////////if (lvl<9) name = GetIonName(Z, A, lvl);
  if (lvl==0 && flb==G4Ions::G4FloatLevelBase::no_Float) name = GetIonName(Z, A, lvl);
  else       name = GetIonName(Z, A, Eex, flb);

  // PDG encoding 
  G4int encoding = GetNucleusEncoding(Z,A,E,lvl);

//G4cout<<"G4IonTable::CreateIon "<<"Z:"<<Z<<" A:"<<A<<" E:"<<E<<" Eex:"<<Eex<<" lvl:"<<lvl<<" name:"<<name<<" code:"<<encoding<<G4endl;
  // PDG mass
  G4double mass =  GetNucleusMass(Z, A)+ Eex;
 
  // PDG charge is set to one of nucleus
  G4double charge =  G4double(Z)*eplus;
 
  // create an ion
  //   spin, parity, isospin values are fixed

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

  //No Anti particle registered
  ion->SetAntiPDGEncoding(0);
  
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4IonTable::CreateIon() : create ion of " << name
	   << "  " << Z << ", " << A
	   << " encoding=" << encoding;
    if (E>0.0) {
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
      for(G4int iCh=0;iCh<nCh;iCh++)
      { decayTable->GetDecayChannel(iCh)->GetDaughter(0); }
    }
  }
#endif
  
  return ion;
}


////////////////////
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, G4int LL, G4double E,
                                            G4Ions::G4FloatLevelBase flb)
{
  if (LL==0) return CreateIon(Z,A,E,flb);
  
  // create hyper nucleus
  G4ParticleDefinition* ion=0;

  // check whether GenericIon has processes
  G4ParticleDefinition* genericIon = 
    G4ParticleTable::GetParticleTable()->GetGenericIon();
  G4ProcessManager* pman=0;
  if (genericIon!=0) pman = genericIon->GetProcessManager();
  if ((genericIon ==0) || (genericIon->GetParticleDefinitionID() < 0) || (pman==0)){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "G4IonTable::CreateIon() : can not create ion of  " 
             << " Z =" << Z << "  A = " << A 
             << "  because GenericIon is not ready !!" <<   G4endl;
    }
#endif
    G4Exception( "G4IonTable::CreateIon()","PART105", JustWarning, 
		 "Can not create ions because GenericIon is not ready");
    return 0;
  }
 
  G4int J=0;
  G4double life = 0.0;
  G4DecayTable* decayTable =0;
  G4bool stable = true;
 
  // excitation energy
  //G4double Eex = G4NuclideTable::Round(E); 
  G4double Eex = E; 
  G4double mass =  GetNucleusMass(Z, A, LL)+ Eex;
  G4int    lvl = 0;
  // lvl is assigned to 9 temporally
  if (Eex>0.0) lvl=9;
   
  // PDG encoding 
  G4int encoding = GetNucleusEncoding(Z,A,LL,E,lvl);

  // PDG charge is set to one of nucleus
  G4double charge =  G4double(Z)*eplus;

  // create an ion
  //   spin, parity, isospin values are fixed
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

  // Release lock for particle table accesses.
  //
 
  G4double mu = 0.0; //  magnetic moment
  ion->SetPDGMagneticMoment(mu);
  static_cast<G4Ions*>(ion)->SetFloatLevelBase(flb);

  //No Anti particle registered
  ion->SetAntiPDGEncoding(0);
  
#ifdef G4VERBOSE
   if (GetVerboseLevel()>1) {
    G4cout << "G4IonTable::CreateIon() : create hyper ion of " << name
	   << "  " << Z << ", " << A << ", " << LL
	   << " encoding=" << encoding;
    if (E>0.0) {
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

////////////////////////////////
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, G4int lvl)
{
  if(lvl == 0) return CreateIon(Z,A,0.0,G4Ions::G4FloatLevelBase::no_Float);
  G4Exception( "G4IonTable::CreateIon()","PART105", JustWarning, 
      "Ion cannot be created by an isomer level. Use excitation energy.");
  return 0;
}


////////////////////
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, G4int LL, G4int lvl)
{
  if (LL==0) return CreateIon(Z,A,lvl);
  if(lvl == 0) return CreateIon(Z,A,0.0,G4Ions::G4FloatLevelBase::no_Float);
  
  if (lvl>0) {
    G4ExceptionDescription ed;
    ed << "Isomer level " << lvl << " is unknown for the isotope (Z="
       << Z << ", A=" << A << ", L=" << LL << "). Null pointer is returned.";
    G4Exception( "G4IonTable::GetIon()","PART106", JustWarning, ed);
    return 0;
  }
  
  return 0;
}

////////////////////
// -- GetIon methods  ------
////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int lvl)
{
  if ( (A<1) || (Z<=0) || (lvl<0) || (A>999) ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A <<  "  Lvl = " << lvl << G4endl;
    }
#endif
    return 0;
  }
  if ( lvl == 0 ) return GetIon(Z,A,0.0);
  
  // Search ions with A, Z, lvl 
  G4ParticleDefinition* ion = FindIon(Z,A,lvl);
  
  // create ion
#ifdef G4MULTITHREADED
  if (ion ==0 ){
    if(G4Threading::IsWorkerThread()){
      G4MUTEXLOCK(&G4IonTable::ionTableMutex);
      ion = FindIonInMaster(Z,A,lvl);
      if(ion != 0) InsertWorker(ion); 
      G4MUTEXUNLOCK(&G4IonTable::ionTableMutex);
    } 
  }
#endif
  if (ion ==0 ){
    G4Exception( "G4IonTable::GetIon()","PART105", JustWarning, 
		 "Ion cannot be created by an isomer level. Use excitation energy.");
    //    G4ExceptionDescription ed;
    //    ed << "Isomer level " << lvl << " is unknown for the isotope (Z="
    //       << Z << ", A=" << A << "). Null pointer is returned.";
    //    G4Exception( "G4IonTable::GetIon()","PART106", JustWarning, ed);
  }
  return ion;  
}


////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int LL, G4int lvl)
{
  if (LL==0) return GetIon(Z,A,lvl);

  if (A < 2 || Z < 0 || Z > A-LL || LL>A || A>999 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A << " L = " << LL 
	     <<"  IsomerLvl = " << lvl << G4endl;
    }
#endif
    return 0;
  } else if( A==2 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : No boud state for " 
             << " Z =" << Z << "  A = " << A << " L = " << LL 
	     <<"  IsomerLvl = " << lvl << G4endl;
    }
#endif
    return 0;
   }

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,LL,lvl);

  // create ion
  if (ion == 0) {
    if (lvl==0) {
#ifdef G4MULTITHREADED
      if(G4Threading::IsWorkerThread()){
        G4MUTEXLOCK(&G4IonTable::ionTableMutex);
        ion = FindIonInMaster(Z,A,LL,lvl);
        if(ion == 0) ion = CreateIon(Z, A, LL, lvl);
        InsertWorker(ion);
        G4MUTEXUNLOCK(&G4IonTable::ionTableMutex);
      } else { 
	ion = CreateIon(Z, A, LL, lvl); 
      }
#else
      ion = CreateIon(Z, A, LL, lvl);
#endif
    } 
  }

//  if(ion == 0)
//  {
//    G4ExceptionDescription ed;
//    ed << "Isomer level " << lvl << " is unknown for the isotope (Z="
//       << Z << ", A=" << A << ", L=" << LL << "). Null pointer is returned.";
//    G4Exception( "G4IonTable::GetIon()","PART106", JustWarning, ed);
//  }
  return ion;  
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4double E, G4int J)
{ return GetIon(Z,A,E,G4Ions::G4FloatLevelBase::no_Float,J); }

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4double E,
                          char flbChar, G4int J)
{ return GetIon(Z,A,E,G4Ions::FloatLevelBase(flbChar),J); }

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4double E,
                          G4Ions::G4FloatLevelBase flb, G4int J)
{
  if ( (A<1) || (Z<=0) || (E<0.0) || (A>999) || (J<0) ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A <<  "  E = " << E/keV << G4endl;
    }
#endif
    return 0;
   }

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,E,flb,J);

  // create ion
#ifdef G4MULTITHREADED
  if(ion == 0)
  {
    if(G4Threading::IsWorkerThread())
    {
      G4MUTEXLOCK(&G4IonTable::ionTableMutex);
      ion = FindIonInMaster(Z,A,E,flb,J);
      if(ion == 0) ion = CreateIon(Z,A,E,flb);
      InsertWorker(ion);
      G4MUTEXUNLOCK(&G4IonTable::ionTableMutex);
    }
    else
    { ion = CreateIon(Z,A,E,flb); }
  }
#else
  if (ion == 0) ion = CreateIon(Z,A,E,flb);
#endif

  return ion;  
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int LL, G4double E, G4int J)
{ return GetIon(Z,A,LL,E,G4Ions::G4FloatLevelBase::no_Float,J); }

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int LL, G4double E,
                           char flbChar, G4int J)
{ return GetIon(Z,A,LL,E,G4Ions::FloatLevelBase(flbChar),J); }

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int LL, G4double E,
                           G4Ions::G4FloatLevelBase flb, G4int J)
{
  if (LL==0) return GetIon(Z,A,E,flb,J);

  if (A < 2 || Z < 0 || Z > A-LL || LL>A || A>999 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A << " L = " << LL 
	     <<"  E = " << E/keV << G4endl;
    }
#endif
    return 0;
  } else if( A==2 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : No boud state for " 
             << " Z =" << Z << "  A = " << A << " L = " << LL 
	     <<  "  E = " << E/keV << G4endl;
    }
#endif
    return 0;
   }

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,LL,E,flb,J);

  // create ion
#ifdef G4MULTITHREADED
  if(ion == 0)
  {
    if(G4Threading::IsWorkerThread())
    {
      G4MUTEXLOCK(&G4IonTable::ionTableMutex);
      ion = FindIonInMaster(Z,A,LL,E,flb,J);
      if(ion == 0) ion = CreateIon(Z,A,LL,E,flb);
      InsertWorker(ion);
      G4MUTEXUNLOCK(&G4IonTable::ionTableMutex);
    }
    else
    { ion = CreateIon(Z,A,LL,E,flb); }
  }
#else
  if(ion == 0) ion = CreateIon(Z,A,LL,E,flb);
#endif

  return ion;  
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int encoding)
{
  G4int Z, A, LL, IsoLvl;
  G4double E;
  if (!GetNucleusByEncoding(encoding,Z,A,LL,E,IsoLvl) ){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : illegal encoding" 
             << " CODE:" << encoding << G4endl;
    }
#endif
    G4Exception( "G4IonTable::GetIon()","PART106",
		 JustWarning, "illegal encoding for an ion");
    return 0;
  }
  //
  return GetIon( Z, A, LL, IsoLvl);
}

/////////////////////
// -- FindIon methods  ------
/////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4double E, G4int J)
{ return FindIon(Z,A,E,G4Ions::G4FloatLevelBase::no_Float,J); }

////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4double E,
                           char flbChar, G4int J)
{ return FindIon(Z,A,E,G4Ions::FloatLevelBase(flbChar),J); }

////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4double E,
                           G4Ions::G4FloatLevelBase flb, G4int J)
{
  if ( (A<1) || (Z<=0) || (J<0) || (E<0.0) || (A>999) ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::FindIon() : illegal atomic number/mass or excitation level " 
             << " Z =" << Z << "  A = " << A <<  "  E = " << E/keV << G4endl;
    }
#endif
    G4Exception( "G4IonTable::FindIon()","PART107",
		 JustWarning, "illegal atomic number/mass");
    return 0;
  }
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion=0;
  G4bool isFound = false;

  // check if light ion
  ion = GetLightIon(Z,A);
  if (ion!=0 && E==0.0) { 
    // light ion 
    isFound = true;
  } else {    
    // -- loop over all particles in Ion table
    G4int encoding=GetNucleusEncoding(Z, A);
    G4IonList::iterator i = fIonList->find(encoding);
    for( ;i != fIonList->end() ; i++) {
      ion = i->second;
      if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
      // excitation level
      G4double anExcitaionEnergy = ((const G4Ions*)(ion))->GetExcitationEnergy();
      if (std::fabs(E - anExcitaionEnergy) < pNuclideTable->GetLevelTolerance() ) {
        if(((const G4Ions*)(ion))->GetFloatLevelBase()==flb)
        {
	  isFound = true;
	  break;
        }
      }
    }
  }

  if ( isFound ){ 
    return const_cast<G4ParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4int LL, G4double E, G4int J)
{ return FindIon(Z,A,LL,E,G4Ions::G4FloatLevelBase::no_Float,J); }

////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4int LL, G4double E,
                         char flbChar, G4int J)
{ return FindIon(Z,A,LL,E,G4Ions::FloatLevelBase(flbChar),J); }

////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4int LL, G4double E,
                         G4Ions::G4FloatLevelBase flb, G4int J)
{
  if (LL==0) return FindIon(Z,A,E,flb,J);
  
  if (A < 2 || Z < 0 || Z > A-LL || LL>A || A>999 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::FindIon() : illegal atomic number/mass or excitation level " 
             << " Z =" << Z << "  A = " << A << " L = " << LL 
	     <<"  E = " << E/keV << G4endl;
    }    
#endif
    G4Exception( "G4IonTable::FindIon()","PART107",
		 JustWarning, "illegal atomic number/mass");
    return 0;
  }
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion=0;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A, LL, 0.0, 0);
  G4IonList::iterator i = fIonList->find(encoding);
  for( ;i != fIonList->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if(  ion->GetQuarkContent(3) != LL) break;
    // excitation level
    G4double anExcitaionEnergy = ((const G4Ions*)(ion))->GetExcitationEnergy();
    if (std::fabs(E - anExcitaionEnergy) < pNuclideTable->GetLevelTolerance() ) {
      if(((const G4Ions*)(ion))->GetFloatLevelBase()==flb)
      {
        isFound = true;
        break;
      }
    }
  }

  if ( isFound ){ 
    return const_cast<G4ParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4int lvl)
{
  if ( (A<1) || (Z<=0) || (lvl<0) || (A>999) ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::FindIon() : illegal atomic number/mass or excitation level " 
             << " Z =" << Z << "  A = " << A <<  "  IsoLvl = " << lvl << G4endl;
    }
#endif
    G4Exception( "G4IonTable::FindIon()","PART107",
		 JustWarning, "illegal atomic number/mass");
    return 0;
  }
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion=0;
  G4bool isFound = false;

  // check if light ion
  ion = GetLightIon(Z,A);
  if (ion!=0 && lvl==0) { 
    // light ion 
    isFound = true;
  } else {    
    // -- loop over all particles in Ion table
    G4int encoding=GetNucleusEncoding(Z, A);
    G4IonList::iterator i = fIonList->find(encoding);
    for( ;i != fIonList->end() ; i++) {
      ion = i->second;
      if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
      // excitation level
      if ( ((const G4Ions*)(ion))->GetIsomerLevel() == lvl) {
	isFound = true;
	break;
      }
    } 
  }

  if ( isFound ){ 
    if(lvl==9)
    {
      G4Exception("G4IonTable::FindIon()","PART5107",JustWarning,
        "Isomer level 9 may be ambiguous.");
    }
    return const_cast<G4ParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4int LL, G4int lvl)
{
  if (LL==0) return FindIon(Z,A,lvl);
  
  if (A < 2 || Z < 0 || Z > A-LL || LL>A || A>999 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::FindIon() : illegal atomic number/mass or excitation level " 
             << " Z =" << Z << "  A = " << A << " L = " << LL 
	     <<"  IsomerLvl = " << lvl << G4endl;
    }    
#endif
    G4Exception( "G4IonTable::FindIon()","PART107",
		 JustWarning, "illegal atomic number/mass");
    return 0;
  }
  // Search ions with A, Z ,E, lvl
  const G4ParticleDefinition* ion=0;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A, LL);
  G4IonList::iterator i = fIonList->find(encoding);
  for( ;i != fIonList->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if ( ion->GetQuarkContent(3) != LL) break;
    // excitation level
    if ( ((const G4Ions*)(ion))->GetIsomerLevel() == lvl) {
      isFound = true;
      break;
    }
  }

  if ( isFound ){ 
    if(lvl==9)
    {
      G4Exception("G4IonTable::FindIon()","PART5107",JustWarning,
        "Isomer level 9 may be ambiguous.");
    }
    return const_cast<G4ParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


/////////////////
G4int G4IonTable::GetNucleusEncoding(G4int Z, G4int A, G4double E, G4int lvl)
{
  // PDG code for Ions
  // Nuclear codes are given as 10-digit numbers +-100ZZZAAAI.
  //For a nucleus consisting of np protons and nn neutrons
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

/////////////////
G4int G4IonTable::GetNucleusEncoding(G4int Z,  G4int A,    G4int LL,
				     G4double E,    G4int lvl)
{
  //  get PDG code for Hyper-Nucleus Ions 
  // Nuclear codes are given as 10-digit numbers +-10LZZZAAAI.
  //For a nucleus consisting of np protons and nn neutrons
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

///////////////
G4bool G4IonTable::GetNucleusByEncoding(G4int encoding,
			    G4int &Z,      G4int &A, 
			    G4double &E,   G4int &lvl)
{
  if (encoding <= 0) return false; // anti particle   

  if (encoding == 2212) {  // proton
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

///////////////
G4bool G4IonTable::GetNucleusByEncoding(G4int encoding,
					G4int &Z,      G4int &A, 
					G4int &LL,   
					G4double &E,   G4int &lvl)
{
  if (encoding <= 0) return false; // anti particle   

 if (encoding == 3122) {  // Lambda
   Z = 1; A = 1; LL = 1;
    E = 0.0; lvl =0; 
    return true;
  } 

  if (encoding % 10 != 0) {
    //!!!not supported for excitation states !!!   
    return false;
  }
  if (encoding < 1000000000) {
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

#include "G4AutoDelete.hh"
/////////////////
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4double E,
                G4Ions::G4FloatLevelBase flb) const 
{
  static G4ThreadLocal G4String *pname = 0;
  if (!pname)  { pname = new G4String(""); G4AutoDelete::Register(pname); }
  G4String &name = *pname;

  static G4ThreadLocal std::ostringstream* os = 0;
  if ( ! os ) {
    os = new std::ostringstream();
    G4AutoDelete::Register(os); 
    os->setf(std::ios::fixed);
    os->precision(3);
  }

  name = GetIonName(Z, A);

  //excited energy
  if ( E>0  || flb!=G4Ions::G4FloatLevelBase::no_Float){
    os->str("");
    std::ostringstream& oo = *os;
    // Excited nucleus
    oo<<'['<<E/keV;
    if(flb!=G4Ions::G4FloatLevelBase::no_Float)
    { oo<<G4Ions::FloatLevelBaseChar(flb); }
    oo<< ']';
    name += os->str();
  }

  return name;
}

/////////////////
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4int LL, G4double E,
                G4Ions::G4FloatLevelBase flb) const 
{
  if (LL==0) return GetIonName(Z, A, E, flb); 
  static G4ThreadLocal G4String *pname = 0;
  if (!pname)  { pname = new G4String(""); G4AutoDelete::Register(pname); }
  G4String &name = *pname;
  name = "";
  for (int i =0; i<LL; i++){
    name +="L";
  }
  name += GetIonName(Z, A, E, flb);
  return name;
}

/////////////////
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4int lvl) const 
{
  static G4ThreadLocal G4String *pname = 0;
  if (!pname)  { pname = new G4String(""); G4AutoDelete::Register(pname); }
  G4String &name = *pname;

  static G4ThreadLocal std::ostringstream* os = 0;
  if ( ! os ) {
    os = new std::ostringstream();
    G4AutoDelete::Register(os);
    os->setf(std::ios::fixed);
  }

  if ( (0< Z) && (Z <=numberOfElements) ) {
    name = elementName[Z-1];
  } else if (Z > numberOfElements) {
    os->str("");
    os->operator<<(Z);
    name = "E" + os->str() + "-";
  } else {
    name = "?";
    return name;
  }
  // Atomic Mass
  os->str("");
  os->operator<<(A);

  if ( lvl>0 ){
    std::ostringstream& oo = *os;
    // isomer level for Excited nucelus 
    oo<<'['<<lvl << ']';
  }
  name += os->str();

  return name;
}

/////////////////
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4int LL, G4int lvl) const 
{
  if (LL==0) return GetIonName(Z, A, lvl); 
  static G4ThreadLocal G4String *pname = 0;
  if (!pname)  { pname = new G4String(""); G4AutoDelete::Register(pname); }
  G4String &name = *pname;
  for (int i =0; i<LL; i++){
    name +="L";
  }
  name += GetIonName(Z, A, lvl);
  return name;
}


/////////////////
G4bool G4IonTable::IsIon(const G4ParticleDefinition* particle)
{
  // return true if the particle is ion

  static const G4String nucleus("nucleus");
  static const G4String proton("proton");

  // neutron is not ion
  if ((particle->GetAtomicMass()>0)   && 
      (particle->GetAtomicNumber()>0) ){
   if (particle->GetBaryonNumber()>0)  return true;
   else return false;
  }

   
  //  particles derived from G4Ions
  if (particle->GetParticleType() == nucleus) return true;

  // proton (Hydrogen nucleus)
  if (particle->GetParticleName() == proton) return true;

  return false;
}

/////////////////
G4bool G4IonTable::IsAntiIon(const G4ParticleDefinition* particle)
{
  // return true if the particle is ion

  static const G4String anti_nucleus("anti_nucleus");
  static const G4String anti_proton("anti_proton");

  // anti_neutron is not ion
  if ((particle->GetAtomicMass()>0)   && 
      (particle->GetAtomicNumber()>0) ){
   if (particle->GetBaryonNumber()<0)  return true;
   else return false;
  }

  //  particles derived from G4Ions
  if (particle->GetParticleType() == anti_nucleus) return true;

  // anti_proton (Anti_Hydrogen nucleus)
  if (particle->GetParticleName() == anti_proton) return true;

  return false;
}

/////////////////
#include <algorithm>

G4bool G4IonTable::IsLightIon(const G4ParticleDefinition* particle) const
{
  static const std::string names[] = { "proton", "alpha", "deuteron",
                           "triton", "He3"};

   // return true if the particle is pre-defined ion
  return std::find(names, names+5, particle->GetParticleName())!=names+5;
} 

G4bool G4IonTable::IsLightAntiIon(const G4ParticleDefinition* particle) const
{
  static const std::string names[] = { "anti_proton", "anti_alpha", "anti_deuteron",
                           "anti_triton", "anti_He3"};

   // return true if the particle is pre-defined ion
  return std::find(names, names+5, particle->GetParticleName())!=names+5;
} 

/////////////////
G4ParticleDefinition* G4IonTable::GetLightIon(G4int Z, G4int A) const
{
  // returns pointer to pre-defined ions
  const G4ParticleDefinition* ion=0;
  if ( (Z<=2) ) {
#ifndef G4MULTITHREADED
      //In sequential use lazy-initialization
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

/////////////////
G4ParticleDefinition* G4IonTable::GetLightAntiIon(G4int Z, G4int A) const
{
  // returns pointer to pre-defined ions 
  const G4ParticleDefinition* ion=0;
  if ( (Z<=2) ) {
#ifndef G4MULTITHREADED
      //In sequential use lazy-initialization
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


/////////////////
// -- GetNucleusMass/GetIonMass ---
/////////////////
G4double  G4IonTable::GetNucleusMass(G4int Z, G4int A, G4int LL, G4int lvl) const
{
  if ( (A<1)  || (Z<0) || (LL<0) || (lvl<0) || (lvl>9) ){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetNucleusMass() : illegal atomic number/mass " 
             << " Z =" << Z << "  A = " << A  
	     << " L = " << LL << " lvl = " << lvl << G4endl;
    }
#endif
    G4Exception( "G4IonTable::GetNucleusMass()","PART107",
		 EventMustBeAborted, "illegal atomic number/mass");
    return -1.0;
  }
  
  G4double mass;
  if (LL == 0) {
    // calculate nucleus mass
    const G4ParticleDefinition* ion=GetLightIon(Z, A);
    
    if (ion!=0) {
      mass = ion->GetPDGMass();
    } else {
      // use G4NucleiProperties::GetNuclearMass
      mass = G4NucleiProperties::GetNuclearMass(A, Z);
    }
    
    // Isomer
    if ( lvl>0 ) {
      // -- loop over all particles in Ion table
      G4int encoding=GetNucleusEncoding(Z, A);
      G4IonList::iterator i = fIonList->find(encoding);
      G4bool isFound =false;
      for( ;i != fIonList->end() ; i++) {
	ion = i->second;
	if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
	// excitation level
	if ( ((const G4Ions*)(ion))->GetIsomerLevel() == lvl) {
	  isFound = true;
	  break;
	}
      }
      if (isFound) {
	// return existing isomer mass
	mass = ion->GetPDGMass();
      } else {
	// Find isomer from IsotopeTable
	const G4IsotopeProperty*  fProperty = FindIsotope(Z, A, lvl);
	if (fProperty !=0 ) mass += fProperty->GetEnergy();
      }
    }

  } else {
    mass = G4HyperNucleiProperties::GetNuclearMass(A, Z, LL);
    }
  return mass;
}

//////////////////
G4double  G4IonTable::GetIsomerMass(G4int Z, G4int A, G4int  lvl) const
{
  return GetNucleusMass(Z,A,0,lvl);
}

//////////////////
G4double  G4IonTable::GetIonMass(G4int Z, G4int A, G4int LL, G4int lvl) const
{
  return GetNucleusMass(Z,A,LL,lvl);
}


/////////////////
// -- Methods for handling conatiner  ---
/////////////////

void G4IonTable::clear()
{
  if (G4ParticleTable::GetParticleTable()->GetReadiness()) {
    G4Exception("G4IonTable::clear()",
		"PART116", JustWarning,
		"No effects because readyToUse is true.");    
    return;
  }

#ifdef G4VERBOSE
    if (GetVerboseLevel()>2) {
      G4cout << "G4IonTable::Clear() : number of Ion regsitered =  "; 
      G4cout << fIonList->size() <<  G4endl;
    }
#endif
  fIonList->clear();
}

void G4IonTable::Insert(const G4ParticleDefinition* particle)
{
  if (!IsIon(particle)) return;
  if (Contains(particle)) return;
 
  G4int Z = particle->GetAtomicNumber();
  G4int A = particle->GetAtomicMass();  
  G4int LL = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, LL); // encoding of the groud state 

  // regsiter the ion with its encoding of the groud state  
  fIonListShadow->insert( std::pair<const G4int, const G4ParticleDefinition*>(encoding, particle) );

}

void G4IonTable::InsertWorker(const G4ParticleDefinition* particle)
{
  if(!particle) return;

  G4int Z = particle->GetAtomicNumber();
  G4int A = particle->GetAtomicMass();  
  G4int LL = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, LL);
  G4bool found = false;
  if (encoding !=0 ) {
    G4IonList::iterator i = fIonList->find(encoding);
    for( ;i != fIonList->end() ; i++) {
      if (particle == i->second ) {
	found  = true;
	break;
      }
    }
  }
  if(found) return;
 
  // regsiter the ion with its encoding of the groud state  
  fIonList->insert( std::pair<const G4int, const G4ParticleDefinition*>(encoding, particle) );

}

/////////////////
void G4IonTable::Remove(const G4ParticleDefinition* particle)
{
  if(!particle) return;
#ifdef G4MULTITHREADED
  if(G4Threading::IsWorkerThread()) {
    G4ExceptionDescription ed;
    ed << "Request of removing " << particle->GetParticleName()
       << " is ignored as it is invoked from a worker thread.";
    G4Exception("G4IonTable::Remove()","PART10117",JustWarning,ed);
    return;
  }
#endif
  if (G4ParticleTable::GetParticleTable()->GetReadiness()) {
    G4StateManager* pStateManager = G4StateManager::GetStateManager();
    G4ApplicationState currentState = pStateManager->GetCurrentState();
    if (currentState != G4State_PreInit) {
      G4String msg = "Request of removing ";
      msg += particle->GetParticleName();  
      msg += " has No effects other than Pre_Init";
      G4Exception("G4IonTable::Remove()",
		  "PART117", JustWarning, msg);
      return;
    } else {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0){
	G4cout << particle->GetParticleName()
	       << " will be removed from the IonTable " << G4endl;
      }
#endif
    }
  }

  if (IsIon(particle)) {
    G4int Z = particle->GetAtomicNumber();
    G4int A = particle->GetAtomicMass();  
    G4int LL = particle->GetQuarkContent(3);  //strangeness 
    G4int encoding=GetNucleusEncoding(Z, A, LL);
    if (encoding !=0 ) {
      G4IonList::iterator i = fIonListShadow->find(encoding);
      for( ;i != fIonListShadow->end() ; i++) {
	if (particle == i->second) {
	  fIonListShadow->erase(i);
	  break;
	}
      }
    }
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "G4IonTable::Remove :" << particle->GetParticleName() 
             << " is not ions" << G4endl; 
    }
#endif
  }
  
}



/////////////////
// -- Dump Information 
/////////////////
void G4IonTable::DumpTable(const G4String &particle_name) const
{
  const G4ParticleDefinition* ion;
  G4IonList::iterator idx;
  for (idx = fIonList->begin(); idx!= fIonList->end(); ++idx) {
    ion = idx->second;
    if (( particle_name == "ALL" ) || (particle_name == "all")){
      ion->DumpTable();
    } else if ( particle_name == ion->GetParticleName() ) {
      ion->DumpTable();
    }
  }
}

/////////////////
const G4String G4IonTable::elementName[] = {
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
              "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", 
              "Cp", "Uut", "Fl","Uup","Lv","Uus","Uuo"
};


/////////////////
G4int G4IonTable::GetVerboseLevel() const
{
  return G4ParticleTable::GetParticleTable()->GetVerboseLevel();
}

/////////////////
void  G4IonTable::AddProcessManager(G4ParticleDefinition* ion)
{
  // check State and do not attach process managaer in event loop
//  G4StateManager* pStateManager = G4StateManager::GetStateManager();
//  G4ApplicationState currentState = pStateManager->GetCurrentState();
//  if (currentState == G4State_EventProc) return;
//  {
//    if (n_error<10)
//    {
//      G4cout << "Defining process manager for " << ion->GetParticleName() << G4endl;
//      G4Exception("G4IonTable::AddProcessManager()", "PART130", JustWarning,
//	"Defining process manager during an event loop is thread unsafe and will be dropped from the next release.");
//      n_error +=1;
//    }
//    return;
//  }

  // check whether GenericIon has processes
  G4ParticleDefinition* genericIon = 
    G4ParticleTable::GetParticleTable()->GetGenericIon();

  G4ProcessManager* pman=0;
  if (genericIon!=0) pman = genericIon->GetProcessManager();
  if ((genericIon ==0) || (genericIon->GetParticleDefinitionID() < 0) || (pman==0)){
    G4cout << "G4IonTable::AddProcessManager() : can not create ion of  " 
           << ion->GetParticleName()
           << "  because GenericIon is not available!!" <<   G4endl;
    G4Exception( "G4IonTable::AddProcessManager()","PART105", FatalException, 
		 "Can not create ions because GenericIon is not available");
    return;
  }
  
////////  ion->SetProcessManager(pman);
  ion->SetParticleDefinitionID(genericIon->GetParticleDefinitionID());
}

#include <vector>     

////////////////////
void  G4IonTable::RegisterIsotopeTable(G4VIsotopeTable* table)
{
  //check duplication
  G4String name = table->GetName();
  for (size_t i = 0; i< fIsotopeTableList->size(); ++i) {
    G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[i];
    if (name == fIsotopeTable->GetName()) return;
  }

  // register 
  fIsotopeTableList->push_back(table);
}

////////////////////
G4VIsotopeTable* G4IonTable::GetIsotopeTable(size_t index) const
{
   G4VIsotopeTable* fIsotopeTable=0;
   if ( index < fIsotopeTableList->size() ) {
     fIsotopeTable = (*fIsotopeTableList)[index];
   }
   return fIsotopeTable;
}


////////////////////
G4IsotopeProperty* G4IonTable::FindIsotope(G4int Z, G4int A, G4double E,
                                     G4Ions::G4FloatLevelBase flb) const
{
  if (fIsotopeTableList ==0) return 0;
  if (fIsotopeTableList->size()==0) return 0;
  
  G4IsotopeProperty* property =0;

  // iterate 
  for (size_t i = 0; i<fIsotopeTableList->size(); ++i) {
    G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[fIsotopeTableList->size()-i-1];
    property = fIsotopeTable->GetIsotope(Z,A,E,flb);
    if(property) break;
  }
  
  return property;
}

////////////////////
G4IsotopeProperty* G4IonTable::FindIsotope(G4int Z, G4int A, G4int lvl) const
{
  if (fIsotopeTableList ==0) return 0;
  if (fIsotopeTableList->size()==0) return 0;
  
  G4IsotopeProperty* property =0;

  // iterate 
  for (size_t i = 0; i<fIsotopeTableList->size(); ++i) {
    G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[fIsotopeTableList->size()-i-1];
    property = fIsotopeTable->GetIsotope(Z,A,lvl);
    if(property) break;
  }
  
  return property;
}


////////////////////
void G4IonTable::CreateAllIon()
{
  PreloadNuclide();
}

////////////////////
void G4IonTable::CreateAllIsomer()
{
  PreloadNuclide();
}

////////////////////
void G4IonTable::PrepareNuclideTable()
{
  if(pNuclideTable==0) pNuclideTable = G4NuclideTable::GetNuclideTable();
}

////////////////////
void G4IonTable::PreloadNuclide()
{
  if ( isIsomerCreated || !G4Threading::IsMultithreadedApplication() ) return;

  pNuclideTable->GenerateNuclide();

  for ( size_t i = 0 ; i != pNuclideTable->entries() ; i++ ) {
     const G4IsotopeProperty*  fProperty = pNuclideTable->GetIsotopeByIndex( i );
     G4int Z  = fProperty->GetAtomicNumber();
     G4int A  = fProperty->GetAtomicMass();
     G4double Eex  = fProperty->GetEnergy();
     GetIon(Z,A,Eex);
  }

  isIsomerCreated = true;
}


////////////////////
G4ParticleDefinition* G4IonTable::GetParticle(G4int index) const
{
  if ( (index >=0) && (index < Entries()) ) {
    G4IonList::iterator idx = fIonList->begin();
    G4int counter = 0;
    while( idx != fIonList->end() ){// Loop checking, 09.08.2015, K.Kurashige
      if ( counter == index ) {
	return const_cast<G4ParticleDefinition*>(idx->second);
      }
      counter++;
      idx++;
    }
  } 
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1){
    G4cout << " G4IonTable::GetParticle"
           << " invalid index (=" << index << ")" 
	   << " entries = " << Entries() << G4endl;
  }
#endif
  return 0;
}

////////////////////
G4bool  G4IonTable::Contains(const G4ParticleDefinition* particle) const
{
  if (!IsIon(particle)) return false;

  G4int Z = particle->GetAtomicNumber();
  G4int A = particle->GetAtomicMass();  
  G4int LL = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, LL);
  G4bool found = false;
  if (encoding !=0 ) {
    G4IonList::iterator i = fIonListShadow->find(encoding);
    for( ;i != fIonListShadow->end() ; i++) {
      if (particle == i->second ) {
	found  = true;
	break;
      }
    }
  }
  return found;
}

////////////////////
G4int G4IonTable::Entries() const
{
  return fIonList->size();
}

////////////////////
G4int G4IonTable::size() const
{
  return fIonList->size();
}


////////////////////
G4ParticleDefinition* G4IonTable::FindIonInMaster(G4int Z, G4int A, G4double E, 
                            G4Ions::G4FloatLevelBase flb, G4int /*J*/)
{
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion=0;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A);
  G4IonList::iterator i = fIonListShadow->find(encoding);
  for( ;i != fIonListShadow->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    // excitation level
    G4double anExcitaionEnergy = ((const G4Ions*)(ion))->GetExcitationEnergy();
    if (std::fabs(E - anExcitaionEnergy) < pNuclideTable->GetLevelTolerance() ) {
      if(((const G4Ions*)(ion))->GetFloatLevelBase()==flb)
      {
        isFound = true;
        break;
      }
    }
  }

  if ( isFound ){ 
    return const_cast<G4ParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
G4ParticleDefinition* G4IonTable::FindIonInMaster(G4int Z, G4int A, G4int LL, G4double E,
                         G4Ions::G4FloatLevelBase flb, G4int J)
{
  if (LL==0) return FindIon(Z,A,E,flb,J);
  
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion=0;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A, LL, 0.0, 0);
  G4IonList::iterator i = fIonListShadow->find(encoding);
  for( ;i != fIonListShadow->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if(  ion->GetQuarkContent(3) != LL) break;
    // excitation level
    G4double anExcitaionEnergy = ((const G4Ions*)(ion))->GetExcitationEnergy();
    if (std::fabs(E - anExcitaionEnergy) < pNuclideTable->GetLevelTolerance() ) {
      if(((const G4Ions*)(ion))->GetFloatLevelBase()==flb)
      {
        isFound = true;
        break;
      }
    }
  }

  if ( isFound ){ 
    return const_cast<G4ParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
G4ParticleDefinition* G4IonTable::FindIonInMaster(G4int Z, G4int A, G4int lvl)
{
  // Search ions with A, Z ,E
  //  !! J is omitted now !!
  const G4ParticleDefinition* ion=0;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A);
  G4IonList::iterator i = fIonListShadow->find(encoding);
  for( ;i != fIonListShadow->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    // excitation level
    if ( ((const G4Ions*)(ion))->GetIsomerLevel() == lvl) {
      isFound = true;
      break;
    } 
  }

  if ( isFound ){ 
    return const_cast<G4ParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
G4ParticleDefinition* G4IonTable::FindIonInMaster(G4int Z, G4int A, G4int LL, G4int lvl)
{
  if (LL==0) return FindIon(Z,A,lvl);
  
  // Search ions with A, Z ,E, lvl
  const G4ParticleDefinition* ion=0;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A, LL);
  G4IonList::iterator i = fIonListShadow->find(encoding);
  for( ;i != fIonListShadow->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if ( ion->GetQuarkContent(3) != LL) break;
    // excitation level
    if ( ((const G4Ions*)(ion))->GetIsomerLevel() == lvl) {
      isFound = true;
      break;
    }
  }

  if ( isFound ){ 
    return const_cast<G4ParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
G4double G4IonTable::GetLifeTime(const G4ParticleDefinition* particle) const
{
  //if(!(particle->IsGeneralIon())) return particle->GetPDGLifeTime();

  //const G4Ions* ion = static_cast<const G4Ions*>(particle);
  //G4int Z = ion->GetAtomicNumber();
  //G4int A = ion->GetAtomicMass();
  //G4double E = ion->GetExcitationEnergy();

  if((particle->IsGeneralIon()) && !pNuclideTable)
  {
   G4Exception("G4IonTable::GetLifeTime()","ParticleIon1001",FatalException,
               "Method is invoked before G4IonTable is initialized.");
   //return 0.;
  } //else {
   //G4IsotopeProperty* isoP = pNuclideTable->GetIsotope(Z,A,E);
    //if(!isoP) return -1001.0;
    //return isoP->GetLifeTime();
  //}
  return particle->GetPDGLifeTime();
}

////////////////////
G4double G4IonTable::GetLifeTime(G4int Z, G4int A, G4double E, char flbChar) const
{ return GetLifeTime(Z,A,E,G4Ions::FloatLevelBase(flbChar)); }

////////////////////
G4double G4IonTable::GetLifeTime(G4int Z, G4int A, G4double E,
             G4Ions::G4FloatLevelBase flb) const
{
  G4double life = -1001.0;
  const G4IsotopeProperty* fProperty = FindIsotope(Z, A, E, flb);
  if( fProperty !=0 ) life = fProperty->GetLifeTime();
  return life;
}


////////////////////
G4ParticleDefinition* G4IonTable::GetMuonicAtom(G4Ions const* base)
{
  if (base==0 || !IsIon(base)){
    G4Exception("G4IonTable::GetMuonicAtom()", "PART987654321",FatalException,"Constructor argument is not a G4Ions");
    return 0;
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
  if(i!=fIonList->end()){
    return const_cast<G4ParticleDefinition*>(i->second);
  }
  // not in threadlocal list; check global list ... 
#ifdef G4MULTITHREADED
  if(G4Threading::IsWorkerThread()){
    G4MUTEXLOCK(&G4IonTable::ionTableMutex);
    i = fIonListShadow->find(encoding);
    auto end = fIonListShadow->end();
    G4MUTEXUNLOCK(&G4IonTable::ionTableMutex);
    if(i!=end){
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
  if(G4Threading::IsWorkerThread()){
    G4MUTEXLOCK(&G4IonTable::ionTableMutex);
    // first, we need to make sure it hasn't been inserted by some
    // other thread
    auto j = fIonListShadow->find(encoding);
    if( j!= fIonListShadow->end() ){
      // oops ... someone else built a copy when we weren't looking;
      // cleanup our instantiation, and take a handle to the one in
      // the global list
      delete muatom;
      muatom = const_cast<G4MuonicAtom*>(static_cast<G4MuonicAtom const*>(j->second));
    } else {
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

////////////////////
G4ParticleDefinition* G4IonTable::GetMuonicAtom(G4int Z, G4int A)
{
  // Need the cast because we need a G4Ions* to pass into the
  // function, but GetIon returns a G4ParticleDefinition* 
  auto base = static_cast<G4Ions const*>(GetIon(Z,A, 0.0));
  return GetMuonicAtom(base);
}




