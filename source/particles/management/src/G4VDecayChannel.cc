// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VDecayChannel.cc,v 1.1 1999-01-07 16:10:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      27 July 1996 H.Kurashige
//      30 May 1997  H.Kurashige
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"

G4VDecayChannel::G4VDecayChannel(const G4String &aName, G4int Verbose)
               :kinematics_name(aName),
		verboseLevel(Verbose),rbranch(0.0),
		numberOfDaughters(0),
		parent_name(NULL),parent(NULL),parent_mass(0.0),
		daughters_name(NULL),daughters(NULL),daughters_mass(NULL),
		particletable(NULL)		
{
  // set pointer to G4ParticleTable (static and singleton object)
  particletable = G4ParticleTable::GetParticleTable();
}

G4VDecayChannel::G4VDecayChannel(const G4String  &aName, 
			       const G4String& theParentName,
			       G4double        theBR,
			       G4int           theNumberOfDaughters,
			       const G4String& theDaughterName1,
			       const G4String& theDaughterName2,
			       const G4String& theDaughterName3,
			       const G4String& theDaughterName4 )
               :kinematics_name(aName),
		verboseLevel(1),rbranch(theBR),
		numberOfDaughters(theNumberOfDaughters),
		parent_name(NULL),parent(NULL),parent_mass(0.0),
		daughters_name(NULL),daughters(NULL),daughters_mass(NULL),
		particletable(NULL)		
{
  // set pointer to G4ParticleTable (static and singleton object)
  particletable = G4ParticleTable::GetParticleTable();

  // parent name
  parent_name = new G4String(theParentName);

  // cleate array
  daughters_name = new G4String*[numberOfDaughters];
  for (G4int index=0;index<numberOfDaughters;index++) daughters_name[index]=NULL;

  // daughters' name
  if (numberOfDaughters>0) daughters_name[0] = new G4String(theDaughterName1);
  if (numberOfDaughters>1) daughters_name[1] = new G4String(theDaughterName2);
  if (numberOfDaughters>2) daughters_name[2] = new G4String(theDaughterName3);
  if (numberOfDaughters>3) daughters_name[3] = new G4String(theDaughterName4);
}



G4VDecayChannel::G4VDecayChannel(const G4VDecayChannel &right)
{
  kinematics_name = right.kinematics_name;
  verboseLevel = right.verboseLevel;
  rbranch = right.rbranch;

  // copy parent name
  parent_name = new G4String(*right.parent_name);
  parent = NULL;
  parent_mass = 0.0; 

  //create array
  numberOfDaughters = right.numberOfDaughters;

  if ( numberOfDaughters >0 ) {
    daughters_name = new G4String*[numberOfDaughters];
    //copy daughters name
    for (G4int index=0; index < numberOfDaughters; index++)
      {
	daughters_name[index] = new G4String(*right.daughters_name[index]);
      }
  }

  //
  daughters_mass = NULL;
  daughters = NULL;

  // particle table
  particletable = G4ParticleTable::GetParticleTable();
}

G4VDecayChannel & G4VDecayChannel::operator=(const G4VDecayChannel &right)
{
  if (this != &right) { 
    kinematics_name = right.kinematics_name;
    verboseLevel = right.verboseLevel;
    rbranch = right.rbranch;

    // copy parent name
    parent_name = new G4String(*right.parent_name);

    // clear daughters_name array
    ClearDaughtersName();

    // recreate array
    numberOfDaughters = right.numberOfDaughters;
    if ( numberOfDaughters >0 ) {
      daughters_name = new G4String*[numberOfDaughters];
      //copy daughters name
      for (G4int index=0; index < numberOfDaughters; index++) {
	  daughters_name[index] = new G4String(*right.daughters_name[index]);
      }
    }
  }

  //
  parent = NULL;
  daughters = NULL;
  parent_mass = 0.0;
  daughters_mass = NULL;

  // particle table
  particletable = G4ParticleTable::GetParticleTable();

  return *this;
}


G4VDecayChannel::~G4VDecayChannel()
{
  if (parent_name != NULL) delete parent_name;
  ClearDaughtersName();
  if (daughters_mass != NULL) delete [] daughters_mass;
} 

void G4VDecayChannel::ClearDaughtersName()
{
  if ( daughters_name != NULL) {
    if (numberOfDaughters>0) {
#ifdef G4VERBOSE
      if (verboseLevel>1) {
	G4cerr << "G4VDecayChannel::ClearDaughtersName ";
	G4cerr << "clear all daughters " << endl;
      }
#endif
      for (G4int index=0; index < numberOfDaughters; index++) { 
	if (daughters_name[index] != NULL) delete daughters_name[index];
      }
    }
    delete [] daughters_name;
    daughters_name = NULL;
  }
  // 
  if (daughters != NULL) delete [] daughters;
  if (daughters_mass != NULL) delete [] daughters_mass;
  daughters = NULL;
  daughters_mass = NULL;

  numberOfDaughters = 0;
}

void G4VDecayChannel::SetNumberOfDaughters(G4int size)
{
  if (size >0) {
    // remove old contents
    ClearDaughtersName();
    // cleate array
    daughters_name = new G4String*[size];
    for (G4int index=0;index<size;index++) daughters_name[index]=NULL;
    numberOfDaughters = size;
  }
}

void G4VDecayChannel::SetDaughter(G4int anIndex, 
				 const G4String &particle_name)
{
  // check numberOfDaughters is positive
  if (numberOfDaughters<=0) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4VDecayChannel::SetDaughter: ";
      G4cerr << "Number of daughters is not defined" << endl;
    }
#endif
    return;
  }

  // check existence of daughters_name array
  if (daughters_name == NULL) {
    // cleate array
    daughters_name = new G4String*[numberOfDaughters];
    for (G4int index=0;index<numberOfDaughters;index++) {
      daughters_name[index]=NULL;
    }
  }

  // check an index    
  if ( (anIndex<0) || (anIndex>=numberOfDaughters) ) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4VDecayChannel::SetDaughter";
      G4cerr << "index out of range " << anIndex << endl;
    }
#endif
  } else {
    // delete the old name if it exists
    if (daughters_name[anIndex]!=NULL) delete daughters_name[anIndex];
    // fill the name
    daughters_name[anIndex] = new G4String(particle_name);
    // refill the array of daughters[] if it exists
    if (daughters != NULL) FillDaughters();
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cerr << "G4VDecayChannel::SetDaughter[" << anIndex <<"] :";
      G4cerr << daughters_name[anIndex] << ":" << *daughters_name[anIndex]<<endl;
    }
#endif
  }
}

void G4VDecayChannel::SetDaughter(G4int anIndex, const G4ParticleDefinition * parent_type)
{
  if (parent_type != NULL) SetDaughter(anIndex, parent_type->GetParticleName());
}

void G4VDecayChannel::FillDaughters()
{
  G4int index;
  
#ifdef G4VERBOSE
  if (verboseLevel>1) G4cerr << "G4VDecayChannel::FillDaughters()" <<endl;
#endif
  if (daughters != NULL) delete [] daughters;

  // parent mass
  if (parent == NULL) FillParent();  
  G4double parentmass = parent->GetPDGMass();

  //
  G4double sumofdaughtermass = 0.0;
  if ((numberOfDaughters <=0) || (daughters_name == NULL) ){
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4VDecayChannel::FillDaughters    ";
      G4cerr << "numberOfDaughters is not defined yet";
    }
#endif
    daughters = NULL;
    G4Exception("G4VDecayChannel::FillDaughters");
  } 
  //create and set the array of pointers to daughter particles
  daughters = new G4ParticleDefinition*[numberOfDaughters];
  if (daughters_mass != NULL) delete [] daughters_mass;
  daughters_mass = new G4double[numberOfDaughters];
  // loop over all daughters
  for (index=0; index < numberOfDaughters;  index++) { 
    if (daughters_name[index] == NULL) {
      // daughter name is not defined
#ifdef G4VERBOSE
      if (verboseLevel>0) {
	G4cerr << "G4VDecayChannel::FillDaughters  ";
	G4cerr << index << "-th daughter is not defined yet" << endl;
      }
#endif
      daughters[index] = NULL;
      G4Exception("G4VDecayChannel::FillDaughters");
    } 
    //search daughter particles in the particle table 
    daughters[index] = particletable->FindParticle(*daughters_name[index]);
    if (daughters[index] == NULL) {
      // can not find the daughter particle
#ifdef G4VERBOSE
      if (verboseLevel>0) {
	G4cerr << "G4VDecayChannel::FillDaughters  ";
	G4cerr << index << ":" << *daughters_name[index];
	G4cerr << " is not defined !!" << endl;
        G4cerr << " The BR of this decay mode is set to zero " << endl;
      }
#endif
      SetBR(0.0);
      // G4Exception("G4VDecayChannel::FillDaughters");
    }
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cerr << index << ":" << *daughters_name[index];
      G4cerr << ":" << daughters[index] << endl;
    }
#endif
    daughters_mass[index] = daughters[index]->GetPDGMass();
    sumofdaughtermass += daughters[index]->GetPDGMass();
  }  // end loop over all daughters

  // check sum of daghter mass
  if (sumofdaughtermass > parentmass) {
    // !!! illegal mass  !!!
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4VDecayChannel::FillDaughters ";
      G4cerr << "    Energy/Momentum conserevation breaks " <<endl;
      G4cerr << "    parent:" << *parent_name;
      G4cerr << " mass:" << parentmass/GeV << "[GeV/c/c]" <<endl;
      for (index=0; index < 3; index++){
        G4cerr << "     daughter " << index << ":" << *daughters_name[index];
        G4cerr << " mass:" << daughters[index]->GetPDGMass()/GeV << "[GeV/c/c]" <<endl;
      }
    }
#endif
    G4Exception("G4VDecayChannel::FillDaughters");
  }
}


void G4VDecayChannel::FillParent()
{
  if (parent_name == NULL) {
    // parent name is not defined
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4VDecayChannel::FillParent   ";
      G4cerr << ": parent name is not defined !!" << endl;
    }
#endif
    parent = NULL;
    G4Exception("G4VDecayChannel::FillParent");
  }
  // search parent particle in the particle table
  parent = particletable->FindParticle(*parent_name);
  if (parent == NULL) {
    // parent particle does not exist
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4VDecayChannel::FillParent   ";
      G4cerr << *parent_name << " does not exist !!" << endl;
    }
#endif
    G4Exception("G4VDecayChannel::FillParent");
  }
  parent_mass = parent->GetPDGMass();
}

void G4VDecayChannel::SetParent(const G4ParticleDefinition * parent_type)
{
  if (parent_type != NULL) SetParent(parent_type->GetParticleName());
}

void G4VDecayChannel::DumpInfo()
{
  G4cout << " BR:  " << rbranch << "  [" << kinematics_name << "]";
  G4cout << "   :  " ;
  for (G4int index=0; index < numberOfDaughters; index++)
  {
    if(daughters_name[index] != NULL) {
      G4cout << " " << *(daughters_name[index]);
    } else {
      G4cout << " not defined ";
    }
  }
  G4cout << endl;
}







