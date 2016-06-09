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

#include "G4ios.hh"
#include <iostream>               
#include <iomanip>               

#include <sstream>


////////////////////
G4IonTable::G4IonTable()
{
  fIonList = new G4IonList();
  fIsotopeTableList = new std::vector<G4VIsotopeTable*>;
}

////////////////////
G4IonTable::~G4IonTable()
{
  // delete IsotopeTable if exists
  if (fIsotopeTableList != 0) {
    for (size_t i = 0; i< fIsotopeTableList->size(); ++i) {
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
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, 
					    G4double E, G4int J)
{
  G4ParticleDefinition* ion=0;

  // check whether the cuurent state is not "PreInit" 
  //  to make sure that GenericIon has processes
  G4ApplicationState currentState = G4StateManager::GetStateManager()->GetCurrentState();
  if (currentState == G4State_PreInit){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cerr << "G4IonTable::CreateIon() : can not create ion of  " 
             << " Z =" << Z << "  A = " << A 
             << "  because the current state is PreInit !!" <<   G4endl;
    }
#endif
    G4Exception( "G4IonTable::CreateIon()","PART105",
		 JustWarning, "Can not create ions in PreInit state");
    return 0;
  }
  
  // get ion name
  G4String name = GetIonName(Z, A, E);
  if ( name(0) == '?') {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4IonTable::CreateIon() : can not create ions " 
             << " Z =" << Z << "  A = " << A <<  G4endl;
    }
#endif
    return 0;
  } 

  G4double life = -1.0;
  G4DecayTable* decayTable =0;
  G4bool stable = true;
  G4double mu = 0.0;

  const G4IsotopeProperty*  fProperty = FindIsotope(Z, A, E, J);
  if (fProperty !=0 ){
    E    = fProperty->GetEnergy();
    J    = fProperty->GetiSpin();
    life = fProperty->GetLifeTime();
    mu   = fProperty->GetMagneticMoment();    
    decayTable = fProperty->GetDecayTable();
  }
  stable = life <= 0.;
  G4double mass =  GetNucleusMass(Z, A)+ E;
  G4double charge =  G4double(Z)*eplus;
 
  G4int encoding = GetNucleusEncoding(Z,A,E,J);

  // create an ion
  //   spin, parity, isospin values are fixed
  //
  ion = new G4Ions(   name,            mass,       0.0*MeV,     charge, 
			 J,              +1,             0,          
			 0,               0,             0,             
		 "nucleus",               0,             A,    encoding,
		    stable,            life,    decayTable,       false,
		  "generic",              0,
		      E                       );
  ion->SetPDGMagneticMoment(mu);

  //No Anti particle registered
  ion->SetAntiPDGEncoding(0);
  
#ifdef G4VERBOSE
   if (GetVerboseLevel()>1) {
    G4cout << "G4IonTable::CreateIon() : create ion of " << name
	   << "  " << Z << ", " << A
	   << " encoding=" << encoding << G4endl;
  } 
#endif
  
  // Add process manager to the ion
  AddProcessManager(name);
 
  return ion;
}


////////////////////
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, G4int L,
					    G4double E, G4int J)
{
  if (L==0) return CreateIon(A,Z,E,J);
  
  // create hyper nucleus
  G4ParticleDefinition* ion=0;

  // check whether the cuurent state is not "PreInit" 
  //  to make sure that GenericIon has processes
  G4ApplicationState currentState = G4StateManager::GetStateManager()->GetCurrentState();
  if (currentState == G4State_PreInit){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cerr << "G4IonTable::CreateIon() : can not create ion of  " 
             << " Z =" << Z << "  A = " << A <<  " L = " <<L 
             << " because the current state is PreInit !!" <<   G4endl;
    }
#endif
    G4Exception( "G4IonTable::CreateIon()","PART105",
		 JustWarning, "Can not create ions in PreInit state");
    return 0;
  }
  
  // get ion name
  G4String name = GetIonName(Z, A, L, E);
  if ( name(L) == '?') {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4IonTable::CreateIon() : can not create ions " 
             << " Z =" << Z << "  A = " << A << "  L = " << L <<  G4endl;
    }
#endif
      return 0;
  } 

  G4double life = -1.0;
  G4DecayTable* decayTable =0;
  G4bool stable = true;
  G4double mu = 0.0;
  G4double mass =  GetNucleusMass(Z, A, L)+ E;
  G4double charge =  G4double(Z)*eplus;
 
  G4int encoding = GetNucleusEncoding(Z,A,L,E,J);

  // create an ion
  //   spin, parity, isospin values are fixed
  //
  ion = new G4Ions(   name,            mass,       0.0*MeV,     charge, 
			 J,              +1,             0,          
			 0,               0,             0,             
		 "nucleus",               0,             A,    encoding,
		    stable,            life,    decayTable,       false,
		  "generic",              0,
		      E                       );
  ion->SetPDGMagneticMoment(mu);

  //No Anti particle registered
  ion->SetAntiPDGEncoding(0);
  
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4IonTable::CreateIon() : create hyper ion of " << name 
	   << " encoding=" << encoding << G4endl;
  } 
#endif
  
  // Add process manager to the ion
  AddProcessManager(name);
      
  return ion;
}

////////////////////
// -- GetIon methods  ------
////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int , G4int )
{
  return GetIon(Z, A);
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int J)
{
  return GetIon( Z, A, 0.0, J);
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int encoding)
{
  G4int Z, A, L, J;
  G4double E;
  if (!GetNucleusByEncoding(encoding,Z,A,L,E,J) ){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4IonTable::GetIon() : illegal encoding" 
             << " CODE:" << encoding << G4endl;
    }
#endif
    G4Exception( "G4IonTable::GetIon()","PART106",
		 JustWarning, "illegal encoding for an ion");
    return 0;
  }
  // Only ground state is supported
  return GetIon( Z, A, L, 0.0, J);
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4double E, G4int J)
{
  if ( (A<1) || (Z<=0) || (J<0) || (E<0.0) || (A>999) ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A <<  "  E = " << E/keV << G4endl;
    }
#endif
    return 0;
   }

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,E,J);

  // create ion
  if (ion == 0) {
    ion = CreateIon(Z, A, E, J);
  }

  return ion;  
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int L, G4double E, G4int J)
{
  if (L==0) return GetIon(Z,A,E,J);

  if (A < 2 || Z < 0 || Z > A-L || L>A || A>999 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4IonTable::GetIon() : illegal atomic number/mass" 
             << " Z =" << Z << "  A = " << A << " L = " << L 
	     <<"  E = " << E/keV << G4endl;
    }
#endif
    return 0;
  } else if( A==2 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4IonTable::GetIon() : No boud state for " 
             << " Z =" << Z << "  A = " << A << " L = " << L 
	     <<  "  E = " << E/keV << G4endl;
    }
#endif
    return 0;
   }

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,L,E,J);

  // create ion
  if (ion == 0) {
    ion = CreateIon(Z, A, L, E, J);
  }

  return ion;  
}

////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4double E, G4int J)
{
  const G4double EnergyTorelance = 0.1 * keV;

  if ( (A<1) || (Z<=0) || (J<0) || (E<0.0) || (A>999) ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4IonTable::FindIon() : illegal atomic number/mass or excitation level " 
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

  // -- loop over all particles in Ion table
  G4int encoding=GetNucleusEncoding(Z, A, 0);
  G4IonList::iterator i = fIonList->find(encoding);
  for( ;i != fIonList->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;

     // excitation level
    G4double anExcitaionEnergy = ((const G4Ions*)(ion))->GetExcitationEnergy();
    if ( ( std::fabs(E - anExcitaionEnergy ) < EnergyTorelance ) ) {
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
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4int L, G4double E, G4int J)
{
  if (L==0) return FindIon(Z,A,E,J);
  
  const G4double EnergyTorelance = 0.1 * keV;

  if (A < 2 || Z < 0 || Z > A-L || L>A || A>999 ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4IonTable::FindIon() : illegal atomic number/mass or excitation level " 
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
  G4int encoding=GetNucleusEncoding(Z, A, L);
  G4IonList::iterator i = fIonList->find(encoding);
  for( ;i != fIonList->end() ; i++) {
    ion = i->second;
    if ( ( ion->GetAtomicNumber() != Z) || (ion->GetAtomicMass()!=A) ) break;
    if(  ion->GetQuarkContent(3) != L) break;

     // excitation level
    G4double anExcitaionEnergy = ((const G4Ions*)(ion))->GetExcitationEnergy();

    if ( ( std::fabs(E - anExcitaionEnergy ) < EnergyTorelance ) ) {
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
G4int G4IonTable::GetNucleusEncoding(G4int Z, G4int A, G4double E, G4int )
{
  // PDG code for Ions
  // Nuclear codes are given as 10-digit numbers +-100ZZZAAAI.
  //For a nucleus consisting of np protons and nn neutrons
  // A = np + nn and Z = np.
  // I gives the isomer level, with I = 0 corresponding 
  // to the ground state and I >0 to excitations
  
  //!!! I = 1 is assigned fo all excitation states !!!   
  const G4double EnergyTorelance = 0.1 * keV;
  if ( Z==1 && A==1 && E< EnergyTorelance ) {
    //proton
    return 2212;
  }
  
  G4int encoding = 1000000000;
  encoding += Z * 10000;
  encoding += A *10;
  if (E>0.0) encoding += 1;
  
  return encoding;
}

/////////////////
G4int G4IonTable::GetNucleusEncoding(G4int Z,  G4int A,    G4int L,
					    G4double E,  G4int )
{
  //  get PDG code for Hyper-Nucleus Ions 
  // Nuclear codes are given as 10-digit numbers +-10LZZZAAAI.
  //For a nucleus consisting of np protons and nn neutrons
  // A = np + nn +nlambda and Z = np.
  // L = nlambda
  // I gives the isomer level, with I = 0 corresponding 
  // to the ground state and I >0 to excitations
  //
  //!!! I = 1 is assigned fo all excitation states in Geant4   

  G4int encoding = 1000000000;
  encoding += L*  10000000;
  encoding += Z * 10000;
  encoding += A *10;
  if (E>0.0) encoding += 1;
  
  return encoding;
}

///////////////
G4bool G4IonTable::GetNucleusByEncoding(G4int encoding,
			    G4int &Z,      G4int &A, 
			    G4double &E,   G4int &J)
{
  if (encoding <= 0) {
    // anti particle   
    return false;
  } 
  if (encoding == 2212) {
    // proton
    Z = 1;
    A = 1;
    E=0.0;
    J=0; 
    return true;
  }

  if (encoding % 10 != 0) {
    //!!!not supported for excitation states !!!   
    return false;
  }

  encoding -= 1000000000;
  Z = encoding/10000;
  encoding -= 10000*Z;
  A = encoding/10;
  
  E=0.0;
  J=0; 
 
  return true;
}
///////////////
G4bool G4IonTable::GetNucleusByEncoding(G4int encoding,
					G4int &Z,      G4int &A, 
					G4int &L,   
					G4double &E,   G4int &J)
{
  if (encoding <= 0) {
    // anti particle   
    return false;
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
  
  E=0.0;
  J=0; 
 
  return true;
}
/////////////////
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4double E) const 
{
  static G4String name;
  name ="";
  if ( (0< Z) && (Z <=numberOfElements) ) {
    name = elementName[Z-1];
  } else if (Z > numberOfElements) {
    std::ostringstream os1;
    os1.setf(std::ios::fixed);
    os1 << Z ;
    name = "E" + os1.str() + "-";
  } else {
    name = "?";
    return name;
  }
  std::ostringstream os;
  os.setf(std::ios::fixed);
  os << A << '[' << std::setprecision(1) << E/keV << ']';
  name += os.str();
  return name;
}

/////////////////
const G4String& G4IonTable::GetIonName(G4int Z, G4int A, G4int L, G4double E) const 
{
  if (L==0) return GetIonName(Z, A, E); 
  static G4String name;
  name ="";
  for (int i =0; i<L; i++){
    name +="L";
  }
  name += GetIonName(Z, A, E);
  return name;
}

/////////////////
G4bool G4IonTable::IsIon(const G4ParticleDefinition* particle)
{
  // return true if the particle is ion

  static G4String nucleus("nucleus");
  static G4String proton("proton");

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

  static G4String anti_nucleus("anti_nucleus");
  static G4String anti_proton("anti_proton");

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
  static G4bool isInitialized = false;
  static const G4ParticleDefinition* p_proton=0;
  static const G4ParticleDefinition* p_deuteron=0;
  static const G4ParticleDefinition* p_triton=0;
  static const G4ParticleDefinition* p_alpha=0;
  static const G4ParticleDefinition* p_He3=0;
  
  if (!isInitialized) {
    p_proton   = G4ParticleTable::GetParticleTable()->FindParticle("proton"); // proton 
    p_deuteron = G4ParticleTable::GetParticleTable()->FindParticle("deuteron"); // deuteron 
    p_triton   = G4ParticleTable::GetParticleTable()->FindParticle("triton"); // tritoon 
    p_alpha    = G4ParticleTable::GetParticleTable()->FindParticle("alpha"); // alpha 
    p_He3      = G4ParticleTable::GetParticleTable()->FindParticle("He3"); // He3 
    isInitialized = true;
  }

  const G4ParticleDefinition* ion=0;
  if ( (Z<=2) ) {
    if ( (Z==1)&&(A==1) ) {
      ion = p_proton;
    } else if ( (Z==1)&&(A==2) ) {
      ion = p_deuteron;
    } else if ( (Z==1)&&(A==3) ) {
      ion = p_triton;
    } else if ( (Z==2)&&(A==4) ) {
      ion = p_alpha;
    } else if ( (Z==2)&&(A==3) ) {
      ion = p_He3;
    }
  }
  return const_cast<G4ParticleDefinition*>(ion);
}

/////////////////
G4ParticleDefinition* G4IonTable::GetLightAntiIon(G4int Z, G4int A) const
{
  // returns pointer to pre-defined ions 
  static G4bool isInitialized = false;
  static const G4ParticleDefinition* p_proton=0;
  static const G4ParticleDefinition* p_deuteron=0;
  static const G4ParticleDefinition* p_triton=0;
  static const G4ParticleDefinition* p_alpha=0;
  static const G4ParticleDefinition* p_He3=0;
  
  if (!isInitialized) {
    p_proton   = G4ParticleTable::GetParticleTable()->FindParticle("anti_proton"); // proton 
    p_deuteron = G4ParticleTable::GetParticleTable()->FindParticle("anti_deuteron"); // deuteron 
    p_triton   = G4ParticleTable::GetParticleTable()->FindParticle("anti_triton"); // tritoon 
    p_alpha    = G4ParticleTable::GetParticleTable()->FindParticle("anti_alpha"); // alpha 
    p_He3      = G4ParticleTable::GetParticleTable()->FindParticle("anti_He3"); // He3 
    isInitialized = true;
  }

  const G4ParticleDefinition* ion=0;
  if ( (Z<=2) ) {
    if ( (Z==1)&&(A==1) ) {
      ion = p_proton;
    } else if ( (Z==1)&&(A==2) ) {
      ion = p_deuteron;
    } else if ( (Z==1)&&(A==3) ) {
      ion = p_triton;
    } else if ( (Z==2)&&(A==4) ) {
      ion = p_alpha;
    } else if ( (Z==2)&&(A==3) ) {
      ion = p_He3;
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
      G4cerr << "G4IonTable::GetNucleusMass() : illegal atomic number/mass " 
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
  G4int encoding=GetNucleusEncoding(Z, A, L);

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
      G4cerr << "G4IonTable::Remove :" << particle->GetParticleName() 
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
              "Uub", "Uut", "Uuq","Uup","Uuh","Uus","Uuo"
};


/////////////////
G4int G4IonTable::GetVerboseLevel() const
{
  return G4ParticleTable::GetParticleTable()->GetVerboseLevel();
}

/////////////////
void  G4IonTable::AddProcessManager(const G4String& name)
{
  // create command string for addProcManager
  std::ostringstream osAdd;
  osAdd << "/run/particle/addProcManager "<< name;
  G4String cmdAdd = osAdd.str();

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
G4IsotopeProperty* G4IonTable::FindIsotope(G4int Z, G4int A, G4double E, G4int )
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
      
#ifdef G4VERBOSE
      if (GetVerboseLevel()>1) {
        G4cout << "G4IonTable::FindIsotope:"; 
        G4cout << " Z: " << Z;
        G4cout << " A: " << A;
        G4cout << " E: " << E;
	G4cout << G4endl; 
	tmp->DumpInfo();
      }
#endif
      if (property !=0) {
	// overwrite spin/magnetic moment/decay table if not defined
	if( property->GetiSpin() ==0) {
	  property->SetiSpin( tmp->GetiSpin() );
	}
	if( property->GetMagneticMoment() <= 0.0) {
	  property->SetMagneticMoment( tmp->GetMagneticMoment() );
	}
	if( property->GetLifeTime() <= 0.0) {
	  property->SetLifeTime( tmp->GetLifeTime() );
	  if (    (property->GetLifeTime() > 0.0)
	       && (property->GetDecayTable() ==0 ) ) {
	    property->SetDecayTable( tmp->GetDecayTable() );
	    tmp->SetDecayTable( 0 );
	  }
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
  G4int    Z; 
  G4int    A;
  G4double E=0.0; 
  G4int    J=0;
  
  for (Z=1; Z<=120; Z++) {
    for (A=Z;A<999 && A<Z*3+10; A++) {
      if (G4NucleiProperties::IsInStableTable(A,Z)){      
	GetIon(Z,A,E,J);
      }
    }
  }
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
    G4cerr << " G4IonTable::GetParticle"
           << " invalid index (=" << index << ")" 
	   << " entries = " << Entries() << G4endl;
  }
#endif
  return 0;
}











