// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VDecayChannel.cc,v 1.7 2000-03-23 16:43:16 hweber Exp $
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
//      23 Mar. 2000 H.Weber      : add GetAngularMomentum
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"

const G4String G4VDecayChannel::noName = " ";

G4VDecayChannel::G4VDecayChannel(const G4String &aName, G4int Verbose)
               :kinematics_name(aName),
		verboseLevel(Verbose),rbranch(0.0),
		numberOfDaughters(0),
		parent_name(0),parent(0),parent_mass(0.0),
		daughters_name(0),daughters(0),daughters_mass(0),
		particletable(0)		
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
		parent_name(0),parent(0),parent_mass(0.0),
		daughters_name(0),daughters(0),daughters_mass(0),
		particletable(0)		
{
  // set pointer to G4ParticleTable (static and singleton object)
  particletable = G4ParticleTable::GetParticleTable();

  // parent name
  parent_name = new G4String(theParentName);

  // cleate array
  daughters_name = new G4String*[numberOfDaughters];
  for (G4int index=0;index<numberOfDaughters;index++) daughters_name[index]=0;

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
  parent = 0;
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
  daughters_mass = 0;
  daughters = 0;

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
  parent = 0;
  daughters = 0;
  parent_mass = 0.0;
  daughters_mass = 0;

  // particle table
  particletable = G4ParticleTable::GetParticleTable();

  return *this;
}


G4VDecayChannel::~G4VDecayChannel()
{
  if (parent_name != 0) delete parent_name;
  ClearDaughtersName();
  if (daughters_mass != 0) delete [] daughters_mass;
} 

void G4VDecayChannel::ClearDaughtersName()
{
  if ( daughters_name != 0) {
    if (numberOfDaughters>0) {
#ifdef G4VERBOSE
      if (verboseLevel>1) {
	G4cout << "G4VDecayChannel::ClearDaughtersName ";
	G4cout << "clear all daughters " << G4endl;
      }
#endif
      for (G4int index=0; index < numberOfDaughters; index++) { 
	if (daughters_name[index] != 0) delete daughters_name[index];
      }
    }
    delete [] daughters_name;
    daughters_name = 0;
  }
  // 
  if (daughters != 0) delete [] daughters;
  if (daughters_mass != 0) delete [] daughters_mass;
  daughters = 0;
  daughters_mass = 0;

  numberOfDaughters = 0;
}

void G4VDecayChannel::SetNumberOfDaughters(G4int size)
{
  if (size >0) {
    // remove old contents
    ClearDaughtersName();
    // cleate array
    daughters_name = new G4String*[size];
    for (G4int index=0;index<size;index++) daughters_name[index]=0;
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
      G4cout << "G4VDecayChannel::SetDaughter: ";
      G4cout << "Number of daughters is not defined" << G4endl;
    }
#endif
    return;
  }

  // check existence of daughters_name array
  if (daughters_name == 0) {
    // cleate array
    daughters_name = new G4String*[numberOfDaughters];
    for (G4int index=0;index<numberOfDaughters;index++) {
      daughters_name[index]=0;
    }
  }

  // check an index    
  if ( (anIndex<0) || (anIndex>=numberOfDaughters) ) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << "G4VDecayChannel::SetDaughter";
      G4cout << "index out of range " << anIndex << G4endl;
    }
#endif
  } else {
    // delete the old name if it exists
    if (daughters_name[anIndex]!=0) delete daughters_name[anIndex];
    // fill the name
    daughters_name[anIndex] = new G4String(particle_name);
    // refill the array of daughters[] if it exists
    if (daughters != 0) FillDaughters();
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << "G4VDecayChannel::SetDaughter[" << anIndex <<"] :";
      G4cout << daughters_name[anIndex] << ":" << *daughters_name[anIndex]<<G4endl;
    }
#endif
  }
}

void G4VDecayChannel::SetDaughter(G4int anIndex, const G4ParticleDefinition * parent_type)
{
  if (parent_type != 0) SetDaughter(anIndex, parent_type->GetParticleName());
}

void G4VDecayChannel::FillDaughters()
{
  G4int index;
  
#ifdef G4VERBOSE
  if (verboseLevel>1) G4cout << "G4VDecayChannel::FillDaughters()" <<G4endl;
#endif
  if (daughters != 0) delete [] daughters;

  // parent mass
  if (parent == 0) FillParent();  
  G4double parentmass = parent->GetPDGMass();

  //
  G4double sumofdaughtermass = 0.0;
  if ((numberOfDaughters <=0) || (daughters_name == 0) ){
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << "G4VDecayChannel::FillDaughters    ";
      G4cout << "numberOfDaughters is not defined yet";
    }
#endif
    daughters = 0;
    G4Exception("G4VDecayChannel::FillDaughters");
  } 
  //create and set the array of pointers to daughter particles
  daughters = new G4ParticleDefinition*[numberOfDaughters];
  if (daughters_mass != 0) delete [] daughters_mass;
  daughters_mass = new G4double[numberOfDaughters];
  // loop over all daughters
  for (index=0; index < numberOfDaughters;  index++) { 
    if (daughters_name[index] == 0) {
      // daughter name is not defined
#ifdef G4VERBOSE
      if (verboseLevel>0) {
	G4cout << "G4VDecayChannel::FillDaughters  ";
	G4cout << index << "-th daughter is not defined yet" << G4endl;
      }
#endif
      daughters[index] = 0;
      G4Exception("G4VDecayChannel::FillDaughters");
    } 
    //search daughter particles in the particle table 
    daughters[index] = particletable->FindParticle(*daughters_name[index]);
    if (daughters[index] == 0) {
      // can not find the daughter particle
#ifdef G4VERBOSE
      if (verboseLevel>0) {
	G4cout << "G4VDecayChannel::FillDaughters  ";
	G4cout << index << ":" << *daughters_name[index];
	G4cout << " is not defined !!" << G4endl;
        G4cout << " The BR of this decay mode is set to zero " << G4endl;
      }
#endif
      SetBR(0.0);
      // G4Exception("G4VDecayChannel::FillDaughters");
    }
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << index << ":" << *daughters_name[index];
      G4cout << ":" << daughters[index] << G4endl;
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
      G4cout << "G4VDecayChannel::FillDaughters ";
      G4cout << "    Energy/Momentum conserevation breaks " <<G4endl;
      G4cout << "    parent:" << *parent_name;
      G4cout << " mass:" << parentmass/GeV << "[GeV/c/c]" <<G4endl;
      for (index=0; index < numberOfDaughters; index++){
        G4cout << "     daughter " << index << ":" << *daughters_name[index];
        G4cout << " mass:" << daughters[index]->GetPDGMass()/GeV << "[GeV/c/c]" <<G4endl;
      }
    }
#endif
    if ( sumofdaughtermass > parentmass+parent->GetPDGWidth() ) {
      G4Exception("G4VDecayChannel::FillDaughters");
    }
  }
}


void G4VDecayChannel::FillParent()
{
  if (parent_name == 0) {
    // parent name is not defined
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << "G4VDecayChannel::FillParent   ";
      G4cout << ": parent name is not defined !!" << G4endl;
    }
#endif
    parent = 0;
    G4Exception("G4VDecayChannel::FillParent");
  }
  // search parent particle in the particle table
  parent = particletable->FindParticle(*parent_name);
  if (parent == 0) {
    // parent particle does not exist
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << "G4VDecayChannel::FillParent   ";
      G4cout << *parent_name << " does not exist !!" << G4endl;
    }
#endif
    G4Exception("G4VDecayChannel::FillParent");
  }
  parent_mass = parent->GetPDGMass();
}

void G4VDecayChannel::SetParent(const G4ParticleDefinition * parent_type)
{
  if (parent_type != 0) SetParent(parent_type->GetParticleName());
}

G4int G4VDecayChannel::GetAngularMomentum()
{
  // determine angular momentum

  // fill pointers to daughter particles if not yet set  
  if (daughters == 0) FillDaughters();

  const G4int PiSpin = parent->GetPDGiSpin();
  const G4int PParity = parent->GetPDGiParity();
  if (2==numberOfDaughters) {     // up to now we can only handle two particle decays
    const G4int D1iSpin  = daughters[0]->GetPDGiSpin();
    const G4int D1Parity = daughters[0]->GetPDGiParity();
    const G4int D2iSpin  = daughters[1]->GetPDGiSpin();
    const G4int D2Parity = daughters[1]->GetPDGiParity();
    const G4int MiniSpin = abs (D1iSpin - D2iSpin);
    const G4int MaxiSpin = D1iSpin + D2iSpin;
    const G4int lMax = (PiSpin+D1iSpin+D2iSpin)/2; // l is allways int
    G4int lMin;
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << "iSpin: " << PiSpin << " -> " << D1iSpin << " + " << D2iSpin << G4endl;
      G4cout << "2*jmin, 2*jmax, lmax " << MiniSpin << " " << MaxiSpin << " " << lMax << G4endl;
    }
#endif
    for (G4int j=MiniSpin; j<=MaxiSpin; j+=2){    // loop over all possible spin couplings
      lMin = abs(PiSpin-j)/2;
#ifdef G4VERBOSE
      if (verboseLevel>1)
	G4cout << "-> checking 2*j=" << j << G4endl;
#endif
      for (G4int l=lMin; l<=lMax; l++) {
#ifdef G4VERBOSE
	if (verboseLevel>1)
	  G4cout << " checking l=" << l << G4endl;
#endif
	if (PParity == D1Parity*D2Parity*pow(-1,l)) {    // check parity for this l
	  return l;
	}
      }
    }
  } else {
    G4Exception ("G4VDecayChannel::GetAngularMomentum: Sorry, can't handle 3 particle decays (up to now)");
  }
  G4Exception ("G4VDecayChannel::GetAngularMomentum: Can't find angular momentum for this decay!");
  return 0;
}

void G4VDecayChannel::DumpInfo()
{
  G4cout << " BR:  " << rbranch << "  [" << kinematics_name << "]";
  G4cout << "   :  " ;
  for (G4int index=0; index < numberOfDaughters; index++)
  {
    if(daughters_name[index] != 0) {
      G4cout << " " << *(daughters_name[index]);
    } else {
      G4cout << " not defined ";
    }
  }
  G4cout << G4endl;
}







