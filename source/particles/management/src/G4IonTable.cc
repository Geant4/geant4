// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IonTable.cc,v 1.11 1999-08-19 08:18:30 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, IT Division, ASD Group
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


#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#include "G4Ions.hh"
#include "G4UImanager.hh"
#include "G4NucleiProperties.hh"

#include "G4ios.hh"

#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif


////////////////////
G4IonTable::G4IonTable()
{
  fIonList = new G4IonList();
}

////////////////////
G4IonTable::~G4IonTable()
{
  G4int idx;

  for (idx=(fIonList->entries()-1); idx >=0 ; idx--) {
    G4ParticleDefinition* particle = (*fIonList)(idx);

    if ( !IsLightIon(particle) ) {
      // delete if not static objects
#ifdef G4VERBOSE
      G4String name;
      if (GetVerboseLevel()>1) {
	G4cout << "G4IonTable:~IonTable() : delete ion of  " ;
        G4cout << particle->GetParticleName() << endl;
      }
#endif
      delete particle;
    }
  }

  // remove all contents in the Ion List 
  fIonList->clear();

  delete fIonList;
}

////////////////////
// -- GetIon methods  ------
////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int J, G4int Q)
{
  return GetIon(Z, A, J);
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4int L)
{

  // Search ions with A, Z 
  G4ParticleDefinition* ion = FindIon(Z,A,L);

  if (ion == 0) {
    // create ion
    G4String name = GetIonName(Z, A, L);
    if ( name(0) == '?') {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0) {
	G4cout << "G4IonTable::GetIon() : can not create ions " << endl;
	G4cout << " Z =" << Z << "  A = " << A <<  endl;
      }
#endif
      return 0;
    } 

    G4double mass =  GetNucleusMass(Z, A);
    G4double charge =  G4double(Z)*eplus;
    G4int J = L; // temporarily spin is set to L/2 (excitaion level)

    // create an ion
    //   spin, parity, isospin values are fixed
    //
    ion = new G4Ions(    name,         mass,       0.0*MeV,     charge, 
			    J,              +1,             0,          
			    0,               0,             0,             
		    "nucleus",       0,             A,           0,
			 true,            -1.0,          0);

#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "G4IonTable::GetIon() : create ion of " << name << endl;
    } 
#endif

    // create command string for addProcManager
    char cmdAdd[60];
    ostrstream osAdd(cmdAdd,60);
    osAdd << "/run/particle/addProcManager "<< name << '\0';
    // set /control/verbose 0
    G4int tempVerboseLevel = G4UImanager::GetUIpointer()->GetVerboseLevel();
    G4UImanager::GetUIpointer()->SetVerboseLevel(0);
    // issue /run/particle/addProcManage
    G4UImanager::GetUIpointer()->ApplyCommand(cmdAdd);
    // retreive  /control/verbose 
    G4UImanager::GetUIpointer()->SetVerboseLevel(tempVerboseLevel);
 
    // Set cut value same as "GenericIon"
    G4ParticleDefinition* genericIon=G4ParticleTable::GetParticleTable()->FindParticle("GenericIon");

    if (genericIon->GetEnergyCuts() != 0) {
      ion->SetCuts( genericIon->GetLengthCuts());
#ifdef G4VERBOSE
      if (GetVerboseLevel()> 1) {
	G4cout << "G4IonTable::GetIon() : cut value =" << genericIon->GetLengthCuts()/mm << "[mm]" <<endl;
      } 
#endif      
      // Build Physics Tables for the ion
      // create command string for buildPhysicsTable
      char cmdBld[60];
      ostrstream osBld(cmdBld,60);
      osBld << "/run/particle/buildPhysicsTable "<< name << '\0';
      // set /control/verbose 0
      tempVerboseLevel = G4UImanager::GetUIpointer()->GetVerboseLevel();
      G4UImanager::GetUIpointer()->SetVerboseLevel(0);
      // issue /run/particle/buildPhysicsTable
      G4UImanager::GetUIpointer()->ApplyCommand(cmdBld);
      // retreive  /control/verbose 
      G4UImanager::GetUIpointer()->SetVerboseLevel(tempVerboseLevel);
    }
  }
  return ion;  
}

////////////////////
G4ParticleDefinition* G4IonTable::GetIon(G4int Z, G4int A, G4double E)
{
   // !!  This method is a dummy now !!
   // !! Implementation will be later !!
   return GetIon(Z, A, 0);  
}

////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4int L)
{
  // Search ions with A, Z 
  G4ParticleDefinition* ion;
  G4bool isFound = false;

  // -- loop over all particles in Ion table
  G4int idx;
  for (idx= 0; idx < fIonList->entries() ; idx++) {
    ion = (*fIonList)(idx);

    // Z = Atomic Number 
    G4int anAtomicNumber = 0;
    // A = baryon number
    G4int anAtomicMass = 0;
    // excitation level
    G4int anExcitaionLevel =0;

    if ( IsLightIon(ion) ) {
      anAtomicNumber = int(ion->GetPDGCharge()/eplus);
      anAtomicMass = ion->GetBaryonNumber();
      anExcitaionLevel =0;
    } else  {
      anAtomicNumber   = ((const G4Ions*)(ion))->GetAtomicNumber();
      anAtomicMass    =  ((const G4Ions*)(ion))->GetAtomicNumber();
      anExcitaionLevel = ((const G4Ions*)(ion))->GetExcitationLevel();
    }

    if ( (A == anAtomicMass) && 
         (Z == anAtomicNumber ) && 
         (L == anExcitaionLevel )            ) {
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

////////////////////
G4ParticleDefinition* G4IonTable::FindIon(G4int Z, G4int A, G4double E)
{
   // !!  This method is a dummy now !!
   // !! Implementation will be later !!
   return FindIon(Z, A, 0);  
}


/////////////////
G4String G4IonTable::GetIonName(G4int Z, G4int A, G4int L) const 
{
  G4String name;
  if ( (0< Z) && (Z <=numberOfElements) ) {
    name = elementName[Z-1];
  } else {
    return "?";
  }
  char val[50];
  ostrstream os(val,50);
  os << A << '[' << L << ']' << '\0';
  name += val;
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
  G4double mass;
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
  // calculate nucleus mass
  G4ParticleDefinition* ion=GetLightIon(A, Z);
  G4double mass;

  if (ion!=0) {
     mass = ion->GetPDGMass();

  } else {

    // This routine returns mass of nuclei (w/o including electron mass) 
    //   mass = Z*proton_mass + (A-Z)*neutron_mass - binding energy

    G4double protonMass = GetProtonMass();
    G4double neutronMass = GetNeutronMass();

    G4double bindingEnergy = G4NucleiProperties::GetNuclearMass(Z, A);
    
     mass = G4double(Z)*protonMass + G4double(A-Z)*neutronMass - bindingEnergy;
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
    fIonList->insert(particle);
  } else {
    //#ifdef G4VERBOSE
    //if (GetVerboseLevel()>0) {
    //  G4cout << "G4IonTable::Insert :" << particle->GetParticleName() ;
    //  G4cout << " is not ions" << endl; 
    //}
    //#endif
  }
}

/////////////////
void G4IonTable::Remove(G4ParticleDefinition* particle)
{
  if (IsIon(particle)) {
    fIonList->remove(particle);
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4IonTable::Remove :" << particle->GetParticleName() ;
      G4cout << " is not ions" << endl; 
    }
#endif
  }
}


/////////////////
// -- Utilities  ---
/////////////////
G4double  G4IonTable::GetProtonMass() const
{
  static G4double protonMass = 0.0;

  // check if proton exits and get the mass
  if (protonMass<=0.0) {
    G4ParticleDefinition* proton = 
          G4ParticleTable::GetParticleTable()->FindParticle("proton");
    if (proton == 0) {
      G4Exception("G4IonTable: G4Proton is not defined !!"); 
    }
    protonMass = proton->GetPDGMass();
  }

  return protonMass;
}

/////////////////
G4double  G4IonTable::GetNeutronMass() const
{
  static G4double neutronMass = 0.0;

  // check if neutron exits and get the mass
  if (neutronMass<=0.0) {
    G4ParticleDefinition* neutron = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
    if (neutron == 0) {
      G4Exception("G4IonTable: G4Neutron is not defined !!"); 
    }
    neutronMass = neutron->GetPDGMass();
  }

  return neutronMass;
}

/////////////////
G4double  G4IonTable::GetElectronMass() const
{
  static G4double electronMass = 0.0;

  // check if electron exits and get the mass
  if (electronMass<=0.0) {
    G4ParticleDefinition* electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    if (electron == 0) {
      G4Exception("G4IonTable: G4Electron is not defined !!"); 
    }
    electronMass = electron->GetPDGMass();
  }

  return electronMass;
}


/////////////////
// -- Dump Information 
/////////////////
void G4IonTable::DumpTable(const G4String &particle_name) const
{
  for (G4int idx= 0; idx < fIonList->entries() ; idx++) {
    if (( particle_name == "ALL" ) || (particle_name == "all")){
      ((*fIonList)(idx))->DumpTable();
    } else if ( particle_name == ((*fIonList)(idx))->GetParticleName() ) {
      ((*fIonList)(idx))->DumpTable();
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
              "La", "Ce", "Pr", "Nd", "Pr", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", 
                   "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", 
  "Fr", "Ra", 
              "Ac", "Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
              "Db", "Jl", "Rf", "Bh", "Hn", "Mt", "Xa"
  
};


/////////////////
G4int G4IonTable::GetVerboseLevel() const
{
  return G4ParticleTable::GetParticleTable()->GetVerboseLevel();
}

































