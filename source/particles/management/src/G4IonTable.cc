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
// $Id: G4IonTable.cc,v 1.41 2006/06/29 19:25:26 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
#include "G4ParticleTable.hh"
#include "G4StateManager.hh"
#include "G4Ions.hh"
#include "G4UImanager.hh"
#include "G4NucleiProperties.hh"

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
  fIsotopeTable = 0;
}

////////////////////
G4IonTable::~G4IonTable()
{
  // delete IsotopeTable if exists
  if (fIsotopeTable != 0) delete fIsotopeTable;
  fIsotopeTable =0;

  if (fIonList ==0) return;

  //  No need to delete here because all particles are dynamic objects
  //   
  // delete ion objects
  //G4ParticleDefinition* particle;
  //G4IonList::reverse_iterator i;
  //for (i = fIonList->rbegin(); i!= fIonList->rend(); ++i) {
  //  particle = *i;
  //
  //if ( !IsLightIon(particle) ) {
  //    delete if not static objects
#ifdef G4VERBOSE
  //    G4String name;
  //    if (GetVerboseLevel()>1) {
  //      G4cout << "G4IonTable:~IonTable() : delete ion of  " ;
  //      G4cout << particle->GetParticleName() << G4endl;
  //    }
#endif
  //    delete particle;
  //  }
  //}

  // remove all contents in the Ion List 
  fIonList->clear();

  delete fIonList;
  fIonList =0;
}


////////////////////
// -- CreateIon method ------
////////////////////
G4ParticleDefinition* G4IonTable::CreateIon(G4int Z, G4int A, G4double E, G4int J)
{
  G4ParticleDefinition* ion=0;

  // check whether the cuurent state is not "PreInit" 
  //  to make sure that GenericIon has processes
  G4ApplicationState currentState = G4StateManager::GetStateManager()->GetCurrentState();
  if (currentState == G4State_PreInit){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "G4IonTable::CreateIon() : can not create ion of  "; 
      G4cout << " Z =" << Z << "  A = " << A <<  G4endl;
      G4cout << " because the current state is PreInit !!" <<   G4endl;
    }
#endif
    G4Exception( "G4IonTable::CreateIon()","Illegal operation",
		 JustWarning, "Can not create ions in PreInit state");
    return 0;
  }
  
  // get ion name
  G4String name = GetIonName(Z, A, E);
  if ( name(0) == '?') {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::CreateIon() : can not create ions " << G4endl;
      G4cout << " Z =" << Z << "  A = " << A <<  G4endl;
    }
#endif
      return 0;
  } 

  G4double life = -1.0;
  G4DecayTable* decayTable =0;
  G4bool stable = true;

  G4IsotopeProperty*  fProperty = FindIsotope(Z, A, E, J);
  if (fProperty !=0 ){
    J =    fProperty->GetiSpin();
    E =    fProperty->GetEnergy();
    life = fProperty->GetLifeTime();
    decayTable = fProperty->GetDecayTable();
  }
  stable = (life <= 0.);
  G4double mass =  GetNucleusMass(Z, A)+ E;
  G4double charge =  G4double(Z)*eplus;
  
  // create an ion
  //   spin, parity, isospin values are fixed
  //
  ion = new G4Ions(   name,            mass,       0.0*MeV,     charge, 
			 J,              +1,             0,          
			 0,               0,             0,             
		 "nucleus",               0,             A,           0,
		    stable,            life,    decayTable,       false,
		  "generic",              0,
		      E                       );

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4IonTable::CreateIon() : create ion of " << name << G4endl;
  } 
#endif
  
  // Add process manager to the ion
  AddProcessManager(name);
  
  // Set cut value same as "GenericIon"
  // SetCuts(ion);
    
  if (fProperty !=0) delete fProperty;
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
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4double E, G4int J)
{
  if ( (A<1) || (Z>numberOfElements) || (Z<=0) || (J<0) || (E<0.0) ) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetIon() : illegal atomic number/mass" << G4endl;
      G4cout << " Z =" << Z << "  A = " << A <<  "  E = " << E/keV << G4endl;
    }
#endif
    G4cerr << "G4IonTable::GetIon called with Z="<<Z<<", A="<<A<<G4endl;
    G4Exception( "G4IonTable::GetIon()","Illegal operation",
		 JustWarning, "illegal atomic number/mass");
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
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4double E, G4int J)
{
  const G4double EnergyTorelance = 0.1 * keV;

  if ( (A<1) || (Z>numberOfElements) || (Z<=0) || (J<0) || (E<0.0)) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::FindIon() : illegal atomic number/mass or excitation level " << G4endl;
      G4cout << " Z =" << Z << "  A = " << A <<  "  E = " << E/keV << G4endl;
    }
#endif
    return 0;
  }
  // Search ions with A, Z ,E
  //  !! J is moitted now !!
  G4ParticleDefinition* ion=0;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4IonList::iterator idx;
  for (idx = fIonList->begin(); idx!= fIonList->end(); ++idx) {
    ion = *idx;

    // Z = Atomic Number 
    G4int     anAtomicNumber   = ion->GetAtomicNumber();
     // A = baryon number
    G4int  anAtomicMass    =  ion->GetAtomicMass();
     // excitation level
    G4double anExcitaionEnergy = ((const G4Ions*)(ion))->GetExcitationEnergy();

    if ( (A == anAtomicMass) && 
         (Z == anAtomicNumber ) && 
	 ( std::abs(E - anExcitaionEnergy ) < EnergyTorelance ) ) {
      isFound = true;
      break;
    }
  }

  if ( isFound ){ 
    return ion;
  } else {
    return 0;
  }
}

/////////////////
G4String G4IonTable::GetIonName(G4int Z, G4int A, G4double E) const 
{
  G4String name;
  if ( (0< Z) && (Z <=numberOfElements) ) {
    name = elementName[Z-1];
  } else {
    return "?";
  }
  std::ostringstream os;
  os.setf(std::ios::fixed);
  os << A << '[' << std::setprecision(1) << E/keV << ']';
  name += os.str();
  return name;
}


/////////////////
G4bool G4IonTable::IsIon(G4ParticleDefinition* particle)
{
  // return true if the particle is ion

  //  particles derived from G4VIon and G4Ions
  G4bool value = (particle->GetParticleType() == "nucleus");

  // proton (Hydrogen nucleus)
  value = value || (particle->GetParticleName() == "proton");

  return value;
}

/////////////////
G4bool G4IonTable::IsLightIon(G4ParticleDefinition* particle) const
{
   // return true if the particle is pre-defined ion
   G4String name = particle->GetParticleName();

   G4bool value  =  (name == "proton");
   value = value || (name == "neutron");
   value = value || (name == "alpha");
   value = value || (name == "deuteron");
   value = value || (name == "triton") ;
   value = value || (name == "He3");
   value = value || (name == "GenericIon") ;
   
   return value;
} 

/////////////////
G4ParticleDefinition* G4IonTable::GetLightIon(G4int Z, G4int A) const
{
  // returns pointer to pre-defined ions 

  G4ParticleDefinition* ion=0;
  if ( (Z<=2) ) {
    if ( (Z==1)&&(A==1) ) {
      ion = G4ParticleTable::GetParticleTable()->FindParticle("proton"); // proton 

    } else if ( (Z==0)&&(A==1) ) {
      ion = G4ParticleTable::GetParticleTable()->FindParticle("neutron"); // neutron 
    } else if ( (Z==1)&&(A==2) ) {
      ion = G4ParticleTable::GetParticleTable()->FindParticle("deuteron"); // deuteron 
    } else if ( (Z==1)&&(A==3) ) {
      ion = G4ParticleTable::GetParticleTable()->FindParticle("triton"); // tritoon 
    } else if ( (Z==2)&&(A==4) ) {
      ion = G4ParticleTable::GetParticleTable()->FindParticle("alpha"); // alpha 

    } else if ( (Z==2)&&(A==3) ) {
      ion = G4ParticleTable::GetParticleTable()->FindParticle("He3"); // He3 
    }
  }
  return ion;
}

/////////////////
// -- GetNucleusMass/GetIonMass ---
/////////////////
G4double  G4IonTable::GetNucleusMass(G4int Z, G4int A) const
{
  if ( (A<1) || (Z>numberOfElements) || (Z<0)) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::GetNucleusMass() : illegal atomic number/mass " << G4endl;
      G4cout << " Z =" << Z << "  A = " << A  << G4endl;
    }
#endif
    G4Exception( "G4IonTable::GetNucleusMass()","Illegal operation",
		 EventMustBeAborted, "illegal atomic number/mass");
    return -1.0;
  }

  // calculate nucleus mass
  G4ParticleDefinition* ion=GetLightIon(Z, A);
  G4double mass;

  if (ion!=0) {
     mass = ion->GetPDGMass();

  } else {

    // use G4NucleiProperties::GetNuclearMass
    mass = G4NucleiProperties::GetNuclearMass(A, Z);

  }

  return mass;
}

//////////////////
G4double  G4IonTable::GetIonMass(G4int Z, G4int A) const
{
   return GetNucleusMass(Z,A);
}


/////////////////
// -- Methods for handling conatiner  ---
/////////////////
void G4IonTable::Insert(G4ParticleDefinition* particle)
{
  if (IsIon(particle)) {
    fIonList->push_back(particle);
  } else {
    //#ifdef G4VERBOSE
    //if (GetVerboseLevel()>0) {
    //  G4cout << "G4IonTable::Insert :" << particle->GetParticleName() ;
    //  G4cout << " is not ions" << G4endl; 
    //}
    //#endif
  }
}

/////////////////
void G4IonTable::Remove(G4ParticleDefinition* particle)
{
  if (IsIon(particle)) {
    G4IonList::iterator idx;
    for (idx = fIonList->begin(); idx!= fIonList->end(); ++idx) {
      if ( particle == *idx) {
        fIonList->erase(idx);
      }
    }
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::Remove :" << particle->GetParticleName() ;
      G4cout << " is not ions" << G4endl; 
    }
#endif
  }

}



/////////////////
// -- Dump Information 
/////////////////
void G4IonTable::DumpTable(const G4String &particle_name) const
{
  G4ParticleDefinition* ion;
  G4IonList::iterator idx;
  for (idx = fIonList->begin(); idx!= fIonList->end(); ++idx) {
    ion = *idx;
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
              "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Xa"
  
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
  fIsotopeTable = table;
}

////////////////////
G4VIsotopeTable* G4IonTable::GetIsotopeTable() const
{
  return fIsotopeTable;
}


////////////////////
G4IsotopeProperty* G4IonTable::FindIsotope(G4int Z, G4int A, G4double E, G4int )
{
  if (fIsotopeTable ==0) return 0;

  // ask IsotopeTable // ask IsotopeTable
  return fIsotopeTable->GetIsotope(Z,A,E);
}














