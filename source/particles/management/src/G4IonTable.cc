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
// $Id$
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
//      Add G4IsomerTable                        5 May. 2013  H.Kurashige

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
#include "G4IsomerTable.hh"

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
////////////////////
G4IonTable::G4IonTable()
  : EnergyUnit (0.1 * keV),
    pIsomerTable(0),
    isIsomerCreated(false)
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
}

// This method is used by each worker thread to copy the content
// from the master thread.
//
void G4IonTable::SlaveG4IonTable()
{
G4Exception("G4ParticleTable::SlaveG4ParticleTable()","G4MT0000",FatalException,"Obsolete");
/*****************************
  fIonList = new G4IonList();
  fIsotopeTableList = new std::vector<G4VIsotopeTable*>;

  G4IonListIterator it;
  for (it = fIonListShadow->begin() ; it != fIonListShadow->end(); it++ )
  {
    fIonList->insert(*it);
  }

  for (size_t i = 0; i < fIsotopeTableListShadow->size(); i++)
  {
    fIsotopeTableList->push_back((*fIsotopeTableListShadow)[i]);
  }
**********************************/
}
void G4IonTable::WorkerG4IonTable()
{
  if( fIonList == 0 )
  { fIonList = new G4IonList(); }
  else
  { fIonList->clear(); }

  G4IonListIterator it;
  for (it = fIonListShadow->begin() ; it != fIonListShadow->end(); it++ )
  {
    G4ParticleDefinition* ion = const_cast<G4ParticleDefinition*>(it->second);
    AddProcessManager(ion); 
    fIonList->insert(*it);
  }

  if( fIsotopeTableList == 0 )
  {
    fIsotopeTableList = new std::vector<G4VIsotopeTable*>; 
    for (size_t i = 0; i < fIsotopeTableListShadow->size(); i++)
    { fIsotopeTableList->push_back((*fIsotopeTableListShadow)[i]); }
  }
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
      delete fIsotopeTable;
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
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, G4double E)
{
  G4ParticleDefinition* ion=0;

  // check whether GenericIon has processes
  G4ParticleDefinition* genericIon = 
    G4ParticleTable::GetParticleTable()->GetGenericIon();
  G4ProcessManager* pman=0;
  if (genericIon!=0) pman = genericIon->GetProcessManager();
  if ((genericIon ==0) || (pman==0)){
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
  
  G4double life = -1.0;
  G4DecayTable* decayTable =0;
  G4bool stable = true;
  G4double mu = 0.0;
  G4double Eex = 0.0;
  G4int    lvl =0;
  G4int    J=0;


  const G4IsotopeProperty*  fProperty = FindIsotope(Z, A, E);
  if (fProperty !=0 ){
    Eex  = fProperty->GetEnergy();
    lvl  = fProperty->GetIsomerLevel();
    J    = fProperty->GetiSpin();
    life = fProperty->GetLifeTime();
    mu   = fProperty->GetMagneticMoment();    
    decayTable = fProperty->GetDecayTable();
    stable = (life <= 0.) || (decayTable ==0);
    lvl = fProperty->GetIsomerLevel();
    if (Eex==0.0) lvl=0;
    if (lvl <0) lvl=9;
  } else {
    // excitation energy
    Eex =  ((int)(E/EnergyUnit+0.5))*EnergyUnit;
    // lvl1 is assigned to 9 temporally    
    if (Eex>0.0) lvl=9;
  }

  // ion name
  G4String name =""; 
  if (lvl<9) name = GetIonName(Z, A, lvl);
  else       name = GetIonName(Z, A, Eex);

  // PDG encoding 
  G4int encoding = GetNucleusEncoding(Z,A,E,lvl);

  // PDG mass
  G4double mass =  GetNucleusMass(Z, A)+ Eex;
 
  // PDG charge is set to one of nucleus
  G4double charge =  G4double(Z)*eplus;
 
  // create an ion
  //   spin, parity, isospin values are fixed


  // Request lock for particle table accesses. Some changes are inside 
  // this critical region.
  //
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4ParticleTable::particleTableMutex);
  G4ParticleTable::lockCount++;
#endif

  ion = new G4Ions(   name,            mass,       0.0*MeV,     charge, 
			 J,              +1,             0,          
			 0,               0,             0,             
		 "nucleus",               0,             A,    encoding,
		    stable,            life,    decayTable,       false,
		 "generic",               0,
		       Eex,             lvl         );

  // Release lock for particle table accesses.
  //
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4ParticleTable::particleTableMutex);
#endif

  ion->SetPDGMagneticMoment(mu);

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
 
  return ion;
}


////////////////////
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, G4int L, G4double E)
{
  if (L==0) return CreateIon(A,Z,E);
  
  // create hyper nucleus
  G4ParticleDefinition* ion=0;

  // check whether GenericIon has processes
  G4ParticleDefinition* genericIon = 
    G4ParticleTable::GetParticleTable()->GetGenericIon();
  G4ProcessManager* pman=0;
  if (genericIon!=0) pman = genericIon->GetProcessManager();
  if ((genericIon ==0) || (pman==0)){
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
  G4double life = -1.0;
  G4DecayTable* decayTable =0;
  G4bool stable = true;
 
  // excitation energy
  G4double Eex =  ((int)(E/EnergyUnit+0.5))*EnergyUnit;
  G4double mass =  GetNucleusMass(Z, A, L)+ Eex;
  G4int    lvl = 0;
  // lvl1 is assigned to 9 temporally
  if (Eex>0.0) lvl=9;
   
  // PDG encoding 
  G4int encoding = GetNucleusEncoding(Z,A,L,E,lvl);

  // PDG charge is set to one of nucleus
  G4double charge =  G4double(Z)*eplus;

  // create an ion
  //   spin, parity, isospin values are fixed
  //
  // get ion name
  G4String name = GetIonName(Z, A, L, Eex);
  
  // Request lock for particle table accesses. Some changes are inside 
  // this critical region.
  //
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4ParticleTable::particleTableMutex);
  G4ParticleTable::lockCount++;
#endif

  ion = new G4Ions(   name,            mass,       0.0*MeV,     charge, 
			 J,              +1,             0,          
			 0,               0,             0,             
		 "nucleus",               0,             A,    encoding,
		    stable,            life,    decayTable,       false,
		  "generic",              0,
		      Eex,              lvl         );

  // Release lock for particle table accesses.
  //
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4ParticleTable::particleTableMutex);
#endif

 
  G4double mu = 0.0; //  magnetic moment
  ion->SetPDGMagneticMoment(mu);

  //No Anti particle registered
  ion->SetAntiPDGEncoding(0);
  
#ifdef G4VERBOSE
   if (GetVerboseLevel()>1) {
    G4cout << "G4IonTable::CreateIon() : create hyper ion of " << name
	   << "  " << Z << ", " << A << ", " << L
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
  G4ParticleDefinition* ion=0;

  // check whether GenericIon has processes
  G4ParticleDefinition* genericIon = 
    G4ParticleTable::GetParticleTable()->GetGenericIon();
  G4ProcessManager* pman=0;
  if (genericIon!=0) pman = genericIon->GetProcessManager();
  if ((genericIon ==0) || (pman==0)){
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
  
  G4double life = -1.0;
  G4DecayTable* decayTable =0;
  G4bool stable = true;
  G4double mu = 0.0;
  G4double Eex=0.0;
  G4int    J=0;

  const G4IsotopeProperty*  fProperty = FindIsotope(Z, A, lvl);
  if (fProperty !=0 ){
    Eex  = fProperty->GetEnergy();
    J    = fProperty->GetiSpin();
    life = fProperty->GetLifeTime();
    mu   = fProperty->GetMagneticMoment();    
    stable = (life <= 0.) || (decayTable ==0);
    lvl = fProperty->GetIsomerLevel();
    if (lvl==0) Eex = 0.0;
    else if (lvl<0) { 
      if (Eex>0.0) lvl = 9;
      else         lvl = 0;
    }
  } else {
    if (lvl>0) return 0;
  }
  // get ion name
  G4String name = GetIonName(Z, A, lvl);
  if (lvl==9) name = GetIonName(Z, A, Eex);

  // PDG encoding 
  G4int encoding = GetNucleusEncoding(Z,A,Eex,lvl);

  // PDG mass
  G4double mass =  GetNucleusMass(Z, A)+ Eex;
 
  // PDG charge is set to one of nucleus
  G4double charge =  G4double(Z)*eplus;
 
  // create an ion
  //   spin, parity, isospin values are fixed
  //

  // Request lock for particle table accesses. Some changes are inside 
  // this critical region.
  //
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4ParticleTable::particleTableMutex);
  G4ParticleTable::lockCount++;
#endif

  ion = new G4Ions(   name,            mass,       0.0*MeV,     charge, 
			 J,              +1,             0,          
			 0,               0,             0,             
		 "nucleus",               0,             A,    encoding,
		    stable,            life,    decayTable,       false,
		 "generic",               0,
		       Eex,             lvl         );

  // Release lock for particle table accesses.
  //
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4ParticleTable::particleTableMutex);
#endif

  ion->SetPDGMagneticMoment(mu);

  //No Anti particle registered
  ion->SetAntiPDGEncoding(0);
  
#ifdef G4VERBOSE
   if (GetVerboseLevel()>1) {
    G4cout << "G4IonTable::CreateIon() : create ion of " << name
	   << "  " << Z << ", " << A
	   << " encoding=" << encoding;
    if (Eex>0.0) {
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


////////////////////
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, G4int L, G4int lvl)
{
  if (L==0) return CreateIon(A,Z,lvl);
  
  // create hyper nucleus
  G4ParticleDefinition* ion=0;

  // check whether GenericIon has processes
  G4ParticleDefinition* genericIon = 
    G4ParticleTable::GetParticleTable()->GetGenericIon();
  G4ProcessManager* pman=0;
  if (genericIon!=0) pman = genericIon->GetProcessManager();
  if ((genericIon ==0) || (pman==0)){
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
 
  if (lvl>0) {
    if (GetVerboseLevel()>1) {
      G4cout << "G4IonTable::CreateIon() : can not create ion of  " 
	     << " Z =" << Z << "  A =" << A  
	     << " L =" << L << " IsomerLvl =" << lvl
	     << "  because excitation energy can not be defined !!" <<   G4endl;
    }
    return 0;
  }
  
  G4int    J = 0; 
  G4double life = -1.0;
  G4DecayTable* decayTable =0;
  G4bool stable = true;
 
  // excitation energy
  G4double Eex =  0.0;
  G4double mass =  GetNucleusMass(Z, A, L)+ Eex;
   
  // PDG charge is set to one of nucleus
  G4double charge =  G4double(Z)*eplus;
  
  // PDG encoding 
  G4int encoding = GetNucleusEncoding(Z,A,L,Eex,lvl);

  // create an ion
  //   spin, parity, isospin values are fixed
  //
  // get ion name
  G4String name = GetIonName(Z, A, L, Eex);

  // Request lock for particle table accesses. Some changes are inside 
  // this critical region.
  //
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4ParticleTable::particleTableMutex);
  G4ParticleTable::lockCount++;
#endif

  ion = new G4Ions(   name,            mass,       0.0*MeV,     charge, 
			 J,              +1,             0,          
			 0,               0,             0,             
		 "nucleus",               0,             A,    encoding,
		    stable,            life,    decayTable,       false,
		  "generic",              0,
		      Eex,              lvl        );

  // Release lock for particle table accesses.
  //
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4ParticleTable::particleTableMutex);
#endif

 
  G4double mu = 0.0; //  magnetic moment
  ion->SetPDGMagneticMoment(mu);

  //No Anti particle registered
  ion->SetAntiPDGEncoding(0);
  
#ifdef G4VERBOSE
   if (GetVerboseLevel()>1) {
    G4cout << "G4IonTable::CreateIon() : create hyper ion of " << name
	   << "  " << Z << ", " << A << ", " << L
	   << " encoding=" << encoding;
    if (Eex>0.0) {
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

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,lvl);

  // create ion
  if (ion == 0) ion = CreateIon(Z, A, lvl);

  return ion;  
}


////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int L, G4int lvl)
{
  if (L==0) return GetIon(Z,A,lvl);

  if (A < 2 || Z < 0 || Z > A-L || L>A || A>999 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A << " L = " << L 
	     <<"  IsomerLvl = " << lvl << G4endl;
    }
#endif
    return 0;
  } else if( A==2 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : No boud state for " 
             << " Z =" << Z << "  A = " << A << " L = " << L 
	     <<"  IsomerLvl = " << lvl << G4endl;
    }
#endif
    return 0;
   }

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,L,lvl);

  // create ion
  if (ion == 0) {
    if (lvl==0) {
      ion = CreateIon(Z, A, L, lvl);
    } else {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
	G4cout << "G4IonTable::GetIon() : can not create ion of  " 
	       << " Z =" << Z << "  A =" << A 
	       << " L= " << L << " IsoLvl=" << lvl 
	       << "  because lvl is non zero  !!" <<   G4endl;
      }
#endif
    }
  }

  return ion;  
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4double E, G4int J)
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
  G4ParticleDefinition* ion = FindIon(Z,A,E,J);

  // create ion
  if (ion == 0) ion = CreateIon(Z, A, E);

  return ion;  
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int L, G4double E, G4int J)
{
  if (L==0) return GetIon(Z,A,E,J);

  if (A < 2 || Z < 0 || Z > A-L || L>A || A>999 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A << " L = " << L 
	     <<"  E = " << E/keV << G4endl;
    }
#endif
    return 0;
  } else if( A==2 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : No boud state for " 
             << " Z =" << Z << "  A = " << A << " L = " << L 
	     <<  "  E = " << E/keV << G4endl;
    }
#endif
    return 0;
   }

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,L,E,J);

  // create ion
  if (ion == 0) ion = CreateIon(Z, A, L, E);

  return ion;  
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int encoding)
{
  G4int Z, A, L, IsoLvl;
  G4double E;
  if (!GetNucleusByEncoding(encoding,Z,A,L,E,IsoLvl) ){
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
  return GetIon( Z, A, L, IsoLvl);
}

////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4double E, G4int J)
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
      if ( std::fabs( E - anExcitaionEnergy )< 10.0*EnergyUnit ) {
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
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4int L, G4double E, G4int J)
{
  if (L==0) return FindIon(Z,A,E,J);
  
  if (A < 2 || Z < 0 || Z > A-L || L>A || A>999 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::FindIon() : illegal atomic number/mass or excitation level " 
             << " Z =" << Z << "  A = " << A << " L = " << L 
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
  G4int encoding=GetNucleusEncoding(Z, A, L, 0.0, 0);
  G4IonList::iterator i = fIonList->find(encoding);
  for( ;i != fIonList->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if(  ion->GetQuarkContent(3) != L) break;
    // excitation level
    G4double anExcitaionEnergy = ((const G4Ions*)(ion))->GetExcitationEnergy();
    if ( std::fabs( E - anExcitaionEnergy )< 10.0*EnergyUnit ) {
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
    return const_cast<G4ParticleDefinition*>(ion);
  } else {
    return 0;
  }
}


////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4int L, G4int lvl)
{
  if (L==0) return FindIon(Z,A,lvl);
  
  if (A < 2 || Z < 0 || Z > A-L || L>A || A>999 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::FindIon() : illegal atomic number/mass or excitation level " 
             << " Z =" << Z << "  A = " << A << " L = " << L 
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
  G4int encoding=GetNucleusEncoding(Z, A, L);
  G4IonList::iterator i = fIonList->find(encoding);
  for( ;i != fIonList->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if ( ion->GetQuarkContent(3) != L) break;
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
  if (lvl>0) encoding +=lvl;       //isomer level
  else if (E>0.0) encoding += 9;  //isomer level
  
  return encoding;
}

/////////////////
G4int G4IonTable::GetNucleusEncoding(G4int Z,  G4int A,    G4int L,
				     G4double E,    G4int lvl)
{
  //  get PDG code for Hyper-Nucleus Ions 
  // Nuclear codes are given as 10-digit numbers +-10LZZZAAAI.
  //For a nucleus consisting of np protons and nn neutrons
  // A = np + nn +nlambda and Z = np.
  // L = nlambda
  // I gives the isomer level, with I = 0 corresponding 
  // to the ground state and I >0 to excitations

  G4int encoding = GetNucleusEncoding(Z, A, E, lvl);
  if (L==0) return encoding; 
  encoding += L*  10000000;
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
					G4int &L,   
					G4double &E,   G4int &lvl)
{
  if (encoding <= 0) return false; // anti particle   

 if (encoding == 3122) {  // Lambda
   Z = 1; A = 1; L = 1;
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
  L = encoding/10000000;
  encoding -= 10000000*L;
  Z = encoding/10000;
  encoding -= 10000*Z;
  A = encoding/10;
  lvl = encoding % 10;
  return true;
}

/////////////////
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4double E) const 
{
  static G4ThreadLocal G4String *pname = 0;
  if (!pname)  { pname = new G4String(""); }
  G4String &name = *pname;

  static G4ThreadLocal std::ostringstream* os = 0;
  if ( ! os ) {
    os = new std::ostringstream();
    os->setf(std::ios::fixed);
    os->precision(1);
  }

  name = GetIonName(Z, A);

  //excited energy
  if ( E>0 ){
    os->str("");
    std::ostringstream& oo = *os;
    // Excited nucelus
    oo<<'['<<E/keV << ']';
    name += os->str();
  }

  return name;
}

/////////////////
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4int L, G4double E) const 
{
  if (L==0) return GetIonName(Z, A, E); 
  static G4ThreadLocal G4String *pname = 0;
  if (!pname)  { pname = new G4String(""); }
  G4String &name = *pname;
  name = "";
  for (int i =0; i<L; i++){
    name +="L";
  }
  name += GetIonName(Z, A, E);
  return name;
}

/////////////////
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4int lvl) const 
{
  static G4ThreadLocal G4String *pname = 0;
  if (!pname)  { pname = new G4String(""); }
  G4String &name = *pname;

  static G4ThreadLocal std::ostringstream* os = 0;
  if ( ! os ) {
    os = new std::ostringstream();
    os->setf(std::ios::fixed);
    os->precision(1);
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
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4int L, G4int lvl) const 
{
  if (L==0) return GetIonName(Z, A, lvl); 
  static G4ThreadLocal G4String *pname = 0;
  if (!pname)  { pname = new G4String(""); }
  G4String &name = *pname;
  for (int i =0; i<L; i++){
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
G4double  G4IonTable::GetNucleusMass(G4int Z, G4int A, G4int L) const
{
  if ( (A<1)  || (Z<0) || (L<0) ){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetNucleusMass() : illegal atomic number/mass " 
             << " Z =" << Z << "  A = " << A  << G4endl;
    }
#endif
    G4Exception( "G4IonTable::GetNucleusMass()","PART107",
		 EventMustBeAborted, "illegal atomic number/mass");
    return -1.0;
  }

  G4double mass;
  if (L == 0) {
    // calculate nucleus mass
    const G4ParticleDefinition* ion=GetLightIon(Z, A);
    
    if (ion!=0) {
      mass = ion->GetPDGMass();
    } else {
      // use G4NucleiProperties::GetNuclearMass
      mass = G4NucleiProperties::GetNuclearMass(A, Z);
    }

  } else {
    mass = G4HyperNucleiProperties::GetNuclearMass(A, Z, L);
  }
  return mass;
}

//////////////////
G4double  G4IonTable::GetIonMass(G4int Z, G4int A, G4int L) const
{
   return GetNucleusMass(Z,A,L);
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
  G4int L = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, L); // encoding of the groud state 

  // regsiter the ion with its encoding of the groud state  
  fIonList->insert( std::pair<const G4int, const G4ParticleDefinition*>(encoding, particle) );

}

/////////////////
void G4IonTable::Remove(const G4ParticleDefinition* particle)
{
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
    G4int L = particle->GetQuarkContent(3);  //strangeness 
    G4int encoding=GetNucleusEncoding(Z, A, L);
    if (encoding !=0 ) {
      G4IonList::iterator i = fIonList->find(encoding);
      for( ;i != fIonList->end() ; i++) {
	if (particle == i->second) {
	  fIonList->erase(i);
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
              "Cn", "Uut", "Fl","Uup","Lv","Uus","Uuo"
};


/////////////////
G4int G4IonTable::GetVerboseLevel() const
{
  return G4ParticleTable::GetParticleTable()->GetVerboseLevel();
}

/////////////////
void  G4IonTable::AddProcessManager(G4ParticleDefinition* ion)
{
  if(ion->GetParticleSubType()!="generic")
  { AddProcessManager(ion->GetParticleName()); return; }
  if(ion->GetParticleName()=="GenericIon")
  { AddProcessManager(ion->GetParticleName()); return; }
  
  // check whether GenericIon has processes
  G4ParticleDefinition* genericIon = 
    G4ParticleTable::GetParticleTable()->GetGenericIon();
  G4ProcessManager* pman=0;
  if (genericIon!=0) pman = genericIon->GetProcessManager();
  if ((genericIon ==0) || (pman==0)){
    G4cout << "G4IonTable::AddProcessManager() : can not create ion of  " 
           << ion->GetParticleName()
           << "  because GenericIon is not ready !!" <<   G4endl;
    G4Exception( "G4IonTable::AddProcessManager()","PART105", FatalException, 
		 "Can not create ions because GenericIon is not ready");
    return;
  }
 
  ion->SetProcessManager(pman);
////  G4cout << "G4IonTable::AddProcessManager() : " << ion->GetParticleName() << G4endl;

}

void  G4IonTable::AddProcessManager(const G4String& name)
{
  // create command string for addProcManager
  std::ostringstream osAdd;
  osAdd << "/run/particle/addProcManager "<< name;
  G4String cmdAdd = osAdd.str();
//G4cout << "G4IonTable::AddProcessManager " << cmdAdd << G4endl;
  // set /control/verbose 0
  G4int tempVerboseLevel = G4UImanager::GetUIpointer()->GetVerboseLevel();
  G4UImanager::GetUIpointer()->SetVerboseLevel(0);

  // issue /run/particle/addProcManage
  G4UImanager::GetUIpointer()->ApplyCommand(cmdAdd);

  // retreive  /control/verbose 
  G4UImanager::GetUIpointer()->SetVerboseLevel(tempVerboseLevel);
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
G4IsotopeProperty* G4IonTable::FindIsotope(G4int Z, G4int A, G4double E)
{
  if (fIsotopeTableList ==0) return 0;
  if (fIsotopeTableList->size()==0) return 0;
  
  // ask IsotopeTable 
  G4IsotopeProperty* property =0;

  // iterate 
  for (size_t i = 0; i< fIsotopeTableList->size(); ++i) {
    G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[i];
    G4IsotopeProperty* tmp = fIsotopeTable->GetIsotope(Z,A,E);
    if ( tmp !=0) {
      G4double Eex  = tmp->GetEnergy(); 
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
        G4cout << "G4IonTable::FindIsotope:"
	       << " Z: " << Z  << " A: " << A
	       << " E: " << E  << G4endl; 
	tmp->DumpInfo();
      }
#endif
      if ( (property !=0) && 
	   ( std::abs(Eex-property->GetEnergy())< 10.0*EnergyUnit)) {
	// overwrite spin/magnetic moment/decay table if not defined
	if( property->GetiSpin() ==0) {
	  property->SetiSpin(tmp->GetiSpin());
	}
	if( property->GetMagneticMoment() <= 0.0) {
	  property->SetMagneticMoment(tmp->GetMagneticMoment());
	}
	if( property->GetLifeTime() <= 0.0) {
	  property->SetLifeTime(tmp->GetLifeTime() );
	  if (    (property->GetLifeTime() > 0.0)
	       && (property->GetDecayTable() ==0 ) ) {
	    property->SetDecayTable( tmp->GetDecayTable() );
	    tmp->SetDecayTable( 0 );
	  }
	}
	if( property->GetIsomerLevel() ==0) {
	  property->SetIsomerLevel(tmp->GetIsomerLevel());
	}
      } else {
	property = tmp;
      }
    }
  }
  
  return property;
}

////////////////////
G4IsotopeProperty* G4IonTable::FindIsotope(G4int Z, G4int A, G4int lvl)
{
  if (fIsotopeTableList ==0) return 0;
  if (fIsotopeTableList->size()==0) return 0;
  
  // ask IsotopeTable 
  G4IsotopeProperty* property =0;

  // iterate 
  for (size_t i = 0; i< fIsotopeTableList->size(); ++i) {
    G4VIsotopeTable* fIsotopeTable= (*fIsotopeTableList)[i];
    G4IsotopeProperty* tmp = fIsotopeTable->GetIsotopeByIsoLvl(Z,A,lvl);
    if ( tmp !=0) {
      G4double Eex  = tmp->GetEnergy();       
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
        G4cout << "G4IonTable::FindIsotope:"
	       << " Z: " << Z  << " A: " << A
	       << " Lvl: " << lvl  << G4endl; 
	tmp->DumpInfo();
      }
#endif
      if (property !=0) {
	// overwrite spin/magnetic moment/decay table if not defined
	if( property->GetiSpin() ==0) {
	  property->SetiSpin(tmp->GetiSpin());
	}
	if( property->GetMagneticMoment() <= 0.0) {
	  property->SetMagneticMoment(tmp->GetMagneticMoment());
	}
	if( property->GetLifeTime() <= 0.0) {
	  property->SetLifeTime(tmp->GetLifeTime() );
	  if (    (property->GetLifeTime() > 0.0)
	       && (property->GetDecayTable() ==0 ) ) {
	    property->SetDecayTable( tmp->GetDecayTable() );
	    tmp->SetDecayTable( 0 );
	  }
	}
	if( (lvl >0) && (property->GetEnergy() ==0.0)) {
	  property->SetEnergy(Eex);
	}
      } else {
	property = tmp;
      }
    }
  }
  
  return property;
}


////////////////////
void G4IonTable::CreateAllIon()
{
  if (pIsomerTable!=0) return;
  
  pIsomerTable = new G4IsomerTable();
  RegisterIsotopeTable(pIsomerTable);

  for (G4int Z=1; Z<=120; Z++) {
    for (G4int A=Z;A<999 && A<Z*3+10; A++) {
      if (G4NucleiProperties::IsInStableTable(A,Z)){      
        GetIon(Z,A);
      }
    }
  }
}

////////////////////
void G4IonTable::CreateAllIsomer()
{
  if (isIsomerCreated) return;
  
  if (pIsomerTable==0) {
    pIsomerTable = new G4IsomerTable();
    RegisterIsotopeTable(pIsomerTable);
  }

  for (G4int Z=1; Z<=120; Z++) {
    for (G4int A=Z;A<999 && A<Z*3+10; A++) {
      if (G4NucleiProperties::IsInStableTable(A,Z)){
	G4int lvl=0;
	GetIon(Z,A,lvl);
	for (lvl=1; lvl<9; lvl++){
	  if( GetIon(Z,A,lvl) ==0) break;;
	}
      }
    }
  }
  isIsomerCreated = true;
}

////////////////////
G4ParticleDefinition* G4IonTable::GetParticle(G4int index) const
{
  if ( (index >=0) && (index < Entries()) ) {
    G4IonList::iterator idx = fIonList->begin();
    G4int counter = 0;
    while( idx != fIonList->end() ){
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
  G4int L = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, L);
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










