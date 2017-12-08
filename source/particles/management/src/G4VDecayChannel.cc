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
// $Id: G4VDecayChannel.cc 105720 2017-08-16 12:38:10Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      27 July 1996 H.Kurashige
//      30 May 1997  H.Kurashige
//      23 Mar. 2000 H.Weber      : add GetAngularMomentum
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4AutoLock.hh"

const G4String G4VDecayChannel::noName = " ";

G4VDecayChannel::G4VDecayChannel()
  :kinematics_name(""),
   rbranch(0.0),
   numberOfDaughters(0),
   parent_name(nullptr), 
   daughters_name(nullptr),
   rangeMass(2.5),
   parent_polarization(),
   particletable(nullptr),
   verboseLevel(1)		
{
  G4MT_parent = nullptr;
  G4MT_daughters = nullptr;
  G4MT_parent_mass = 0.0;
  G4MT_daughters_mass = nullptr;
  G4MT_daughters_width = nullptr;

  // set pointer to G4ParticleTable (static and singleton object)
  particletable = G4ParticleTable::GetParticleTable();
}

G4VDecayChannel::G4VDecayChannel(const G4String &aName, G4int Verbose)
  :kinematics_name(aName),
   rbranch(0.0),
   numberOfDaughters(0),
   parent_name(nullptr), 
   daughters_name(nullptr),
   rangeMass(2.5),
   parent_polarization(),
   particletable(nullptr),
   verboseLevel(Verbose),
   daughtersMutex(G4MUTEX_INITIALIZER),
   parentMutex(G4MUTEX_INITIALIZER)
{
  G4MT_parent = nullptr;
  G4MT_daughters = nullptr;
  G4MT_parent_mass = 0.0;
  G4MT_daughters_mass = nullptr;
  G4MT_daughters_width = nullptr;

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
		rbranch(theBR),
		numberOfDaughters(theNumberOfDaughters),
		parent_name(nullptr), 
		daughters_name(nullptr),
                rangeMass(1.0),
                parent_polarization(),
		particletable(0),
		verboseLevel(1),
		daughtersMutex(G4MUTEX_INITIALIZER),
		parentMutex(G4MUTEX_INITIALIZER)
{
  G4MT_parent = nullptr;
  G4MT_daughters = nullptr;
  G4MT_parent_mass = 0.0;
  G4MT_daughters_mass = nullptr;
  G4MT_daughters_width = nullptr;

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

  if      (rbranch <0.  ) rbranch = 0.0;
  else if (rbranch >1.0 ) rbranch = 1.0;
}

G4VDecayChannel::G4VDecayChannel(const G4VDecayChannel &right)
{
  kinematics_name = right.kinematics_name;
  verboseLevel = right.verboseLevel;
  rbranch = right.rbranch;
  rangeMass =  right.rangeMass;

  // copy parent name
  parent_name = new G4String(*right.parent_name);
  G4MT_parent = nullptr;
  G4MT_parent_mass = 0.0; 

  //create array
  numberOfDaughters = right.numberOfDaughters;

  daughters_name =nullptr;
  if ( numberOfDaughters >0 ) {
    daughters_name = new G4String*[numberOfDaughters];
    //copy daughters name
    for (G4int index=0; index < numberOfDaughters; index++){
      daughters_name[index] = new G4String(*right.daughters_name[index]);
    }
  }

  //
  G4MT_daughters_mass = nullptr;
  G4MT_daughters = nullptr;
  G4MT_daughters_width = nullptr;

  // particle table
  particletable = G4ParticleTable::GetParticleTable();

  parent_polarization = right.parent_polarization;

  G4MUTEXINIT(daughtersMutex);
  G4MUTEXINIT(parentMutex);
}

G4VDecayChannel & G4VDecayChannel::operator=(const G4VDecayChannel &right)
{
  if (this != &right) { 
    kinematics_name = right.kinematics_name;
    verboseLevel = right.verboseLevel;
    rbranch = right.rbranch;
    rangeMass =  right.rangeMass;
    parent_polarization = right.parent_polarization;
    // copy parent name
    parent_name = new G4String(*right.parent_name);

    // clear daughters_name array
    ClearDaughtersName();

    // recreate array
    numberOfDaughters = right.numberOfDaughters;
    if ( numberOfDaughters >0 ) {
      if (daughters_name != nullptr) ClearDaughtersName();
      daughters_name = new G4String*[numberOfDaughters];
      //copy daughters name
      for (G4int index=0; index < numberOfDaughters; index++) {
	  daughters_name[index] = new G4String(*right.daughters_name[index]);
      }
    }
  }

  //
  G4MT_parent = nullptr;
  G4MT_daughters = nullptr;
  G4MT_parent_mass = 0.0;
  G4MT_daughters_mass = nullptr;
  G4MT_daughters_width = nullptr;

  // particle table
  particletable = G4ParticleTable::GetParticleTable();

  G4MUTEXINIT(daughtersMutex);
  G4MUTEXINIT(parentMutex);

  return *this;
}

G4VDecayChannel::~G4VDecayChannel()
{
  ClearDaughtersName();
  if (parent_name != nullptr) delete parent_name;
  parent_name = nullptr;
  if (G4MT_daughters_mass != nullptr) delete [] G4MT_daughters_mass;
  G4MT_daughters_mass =nullptr;
  if (G4MT_daughters_width != nullptr) delete [] G4MT_daughters_width;
  G4MT_daughters_width = nullptr;
  G4MUTEXDESTROY(daughtersMutex);
  G4MUTEXDESTROY(parentMutex);
} 

void G4VDecayChannel::ClearDaughtersName()
{
  G4AutoLock l(&daughtersMutex);
  if ( daughters_name != nullptr) {
    if (numberOfDaughters>0) {
#ifdef G4VERBOSE
      if (verboseLevel>1) {
	G4cout << "G4VDecayChannel::ClearDaughtersName "
	       << " for " << *parent_name << G4endl;
      }
#endif
      for (G4int index=0; index < numberOfDaughters; index++) { 
	if (daughters_name[index] != nullptr) delete daughters_name[index];
      }
    }
    delete [] daughters_name;
    daughters_name = nullptr;
  }
  // 
  if (G4MT_daughters != nullptr) delete [] G4MT_daughters;
  if (G4MT_daughters_mass != nullptr) delete [] G4MT_daughters_mass;
  if (G4MT_daughters_width != nullptr) delete [] G4MT_daughters_width;
  G4MT_daughters_width = nullptr;
  G4MT_daughters = nullptr;
  G4MT_daughters_mass = nullptr;

  numberOfDaughters = 0;
}

void G4VDecayChannel::SetNumberOfDaughters(G4int size)
{
  if (size >0) {
    // remove old contents
    ClearDaughtersName();
    // cleate array
    daughters_name = new G4String*[size];
    for (G4int index=0;index<size;index++) daughters_name[index]=nullptr;
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
      G4cout << "G4VDecayChannel::SetDaughter: "
             << "Number of daughters is not defined" << G4endl;
    }
#endif
    return;
  }

  //ANDREA:-> Feb 25 2016
  // An analysis of this code, shows that this method is called
  // only in the constructor of derived classes.
  // The general idea of this method is probably to support
  // the possibility to re-define daughters on the fly, however
  // this design is extremely problematic for MT mode, we thus
  // require (as practically happens) that the method is called only
  // at construction, i.e. when G4MT_daugheters == 0
  // moreover this method can be called only after SetNumberOfDaugthers
  // has been called (see previous if), in such a case daughters_name != 0
  if ( daughters_name == nullptr ) {
    G4Exception("G4VDecayChannel::SetDaughter","PART112",FatalException,
		"Trying to add a daughter without specifying number of secondaries, useSetNumberOfDaughters first");
    return;
  }
  if ( G4MT_daughters != nullptr ) {
    G4Exception("G4VDecayChannel::SetDaughter","PART111",FatalException,
		"Trying to modify a daughter of a decay channel, but decay channel already has daughters.");
    return;
  }
  //<-:ANDREA

  // check an index    
  if ( (anIndex<0) || (anIndex>=numberOfDaughters) ) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << "G4VDecayChannel::SetDaughter"
             << "index out of range " << anIndex << G4endl;
    }
#endif
  } else {
    // fill the name
    daughters_name[anIndex] = new G4String(particle_name);
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
  if (parent_type != nullptr) SetDaughter(anIndex, parent_type->GetParticleName());
}

void G4VDecayChannel::FillDaughters()
{
  G4AutoLock lock(&daughtersMutex);
  //Double check, check again if another thread has already filled this, in
  //case do not need to do anything
  if ( G4MT_daughters != nullptr ) return;

  G4int index;
  
#ifdef G4VERBOSE
  if (verboseLevel>1) G4cout << "G4VDecayChannel::FillDaughters()" <<G4endl;
#endif
  if (G4MT_daughters != nullptr) {
    delete [] G4MT_daughters;
    G4MT_daughters = nullptr;
  }

  // parent mass
  CheckAndFillParent();
  G4double parentmass = G4MT_parent->GetPDGMass();

  //
  G4double sumofdaughtermass = 0.0;
  G4double sumofdaughterwidthsq = 0.0;

  if ((numberOfDaughters <=0) || (daughters_name == nullptr) ){
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << "G4VDecayChannel::FillDaughters   "
             << "[ " << G4MT_parent->GetParticleName() << " ]"
             << "numberOfDaughters is not defined yet";
    }
#endif
    G4MT_daughters = nullptr;
    G4Exception("G4VDecayChannel::FillDaughters",
		"PART011", FatalException,
		"Can not fill daughters: numberOfDaughters is not defined yet");    
  } 

  //create and set the array of pointers to daughter particles
  G4MT_daughters = new G4ParticleDefinition*[numberOfDaughters];
  if (G4MT_daughters_mass != nullptr) delete [] G4MT_daughters_mass;
  if (G4MT_daughters_width != nullptr) delete [] G4MT_daughters_width;
  G4MT_daughters_mass = new G4double[numberOfDaughters];
  G4MT_daughters_width = new G4double[numberOfDaughters];
  // loop over all daughters
  for (index=0; index < numberOfDaughters;  index++) { 
    if (daughters_name[index] == nullptr) {
      // daughter name is not defined
#ifdef G4VERBOSE
      if (verboseLevel>0) {
	G4cout << "G4VDecayChannel::FillDaughters  "
	       << "[ " << G4MT_parent->GetParticleName() << " ]"
	       << index << "-th daughter is not defined yet" << G4endl;
      }
#endif
      G4MT_daughters[index] = nullptr;
      G4Exception("G4VDecayChannel::FillDaughters",
		  "PART011", FatalException,
		  "Can not fill daughters: name of a daughter is not defined yet");    
    } 
    //search daughter particles in the particle table 
    G4MT_daughters[index] = particletable->FindParticle(*daughters_name[index]);
    if (G4MT_daughters[index] == nullptr ) {
      // can not find the daughter particle
#ifdef G4VERBOSE
      if (verboseLevel>0) {
	G4cout << "G4VDecayChannel::FillDaughters  "
	        << "[ " << G4MT_parent->GetParticleName() << " ]"
               << index << ":" << *daughters_name[index]
	       << " is not defined !!" << G4endl;
        G4cout << " The BR of this decay mode is set to zero " << G4endl;
      }
#endif
      SetBR(0.0);
      return;
    }
#ifdef G4VERBOSE
    if (verboseLevel>1) {
      G4cout << index << ":" << *daughters_name[index];
      G4cout << ":" << G4MT_daughters[index] << G4endl;
    }
#endif
    G4MT_daughters_mass[index] = G4MT_daughters[index]->GetPDGMass();
    G4double d_width = G4MT_daughters[index]->GetPDGWidth();
    G4MT_daughters_width[index] = d_width;
    sumofdaughtermass += G4MT_daughters[index]->GetPDGMass();
    sumofdaughterwidthsq += d_width*d_width;
  }  // end loop over all daughters

  // check sum of daghter mass
  G4double widthMass = std::sqrt(G4MT_parent->GetPDGWidth()*G4MT_parent->GetPDGWidth()+sumofdaughterwidthsq);
  if ( (G4MT_parent->GetParticleType() != "nucleus") &&
       (numberOfDaughters !=1) &&
       (sumofdaughtermass > parentmass + rangeMass*widthMass) ){
   // !!! illegal mass  !!!
#ifdef G4VERBOSE
   if (GetVerboseLevel()>0) {
     G4cout << "G4VDecayChannel::FillDaughters "
            << "[ " << G4MT_parent->GetParticleName() << " ]"
            << "    Energy/Momentum conserevation breaks " <<G4endl;
     if (GetVerboseLevel()>1) {
       G4cout << "    parent:" << *parent_name
              << " mass:" << parentmass/GeV << "[GeV/c/c]" <<G4endl;
       for (index=0; index < numberOfDaughters; index++){
	 G4cout << "     daughter " << index << ":" << *daughters_name[index]
	        << " mass:" << G4MT_daughters[index]->GetPDGMass()/GeV
  	        << "[GeV/c/c]" <<G4endl;
       }
     }
   }
#endif
 }
}


void G4VDecayChannel::FillParent()
{
  G4AutoLock lock(&parentMutex);
  //Double check, check again if another thread has already filled this, in
  //case do not need to do anything
  if ( G4MT_parent != nullptr ) return;

  if (parent_name == nullptr) {
    // parent name is not defined
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << "G4VDecayChannel::FillParent   "
             << ": parent name is not defined !!" << G4endl;
    }
#endif
    G4MT_parent = nullptr;
    G4Exception("G4VDecayChannel::FillParent()",
		"PART012", FatalException,
		"Can not fill parent: parent name is not defined yet");    
    return;
  }
  // search parent particle in the particle table
  G4MT_parent = particletable->FindParticle(*parent_name);
  if (G4MT_parent == nullptr) {
    // parent particle does not exist
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cout << "G4VDecayChannel::FillParent   "
             << *parent_name << " does not exist !!" << G4endl;
    }
#endif
    G4Exception("G4VDecayChannel::FillParent()",
		"PART012", FatalException,
		"Can not fill parent: parent does not exist");
    return;
  }
  G4MT_parent_mass = G4MT_parent->GetPDGMass();
}

void G4VDecayChannel::SetParent(const G4ParticleDefinition * parent_type)
{
  if (parent_type != nullptr) SetParent(parent_type->GetParticleName());
}

G4int G4VDecayChannel::GetAngularMomentum()
{
  // determine angular momentum

  // fill pointers to daughter particles if not yet set  
  CheckAndFillDaughters();

  const G4int PiSpin = G4MT_parent->GetPDGiSpin();
  const G4int PParity = G4MT_parent->GetPDGiParity();
  if (2==numberOfDaughters) {     // up to now we can only handle two particle decays
    const G4int D1iSpin  = G4MT_daughters[0]->GetPDGiSpin();
    const G4int D1Parity = G4MT_daughters[0]->GetPDGiParity();
    const G4int D2iSpin  = G4MT_daughters[1]->GetPDGiSpin();
    const G4int D2Parity = G4MT_daughters[1]->GetPDGiParity();
    const G4int MiniSpin = std::abs (D1iSpin - D2iSpin);
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
      lMin = std::abs(PiSpin-j)/2;
#ifdef G4VERBOSE 
      if (verboseLevel>1)
	G4cout << "-> checking 2*j=" << j << G4endl;
#endif
      for (G4int l=lMin; l<=lMax; l++) {
#ifdef G4VERBOSE
	if (verboseLevel>1)
	  G4cout << " checking l=" << l << G4endl;
#endif
        if (l%2==0) {
	  if (PParity == D1Parity*D2Parity) {    // check parity for this l
	    return l;
          } 
	} else {
	  if (PParity == -1*D1Parity*D2Parity) {    // check parity for this l
            return l;
          }
        }
      }
    }
  } else {
    G4Exception("G4VDecayChannel::GetAngularMomentum",
		"PART111", JustWarning,
		"Sorry, can't handle 3 particle decays (up to now)");
    return 0;
  }
  G4Exception ("G4VDecayChannel::GetAngularMomentum",
		"PART111", JustWarning,
		"Can't find angular momentum for this decay");
  return 0;
}

void G4VDecayChannel::DumpInfo()
{
  G4cout << " BR:  " << rbranch << "  [" << kinematics_name << "]";
  G4cout << "   :  " ;
  for (G4int index=0; index < numberOfDaughters; index++){
    if(daughters_name[index] != nullptr) {
      G4cout << " " << *(daughters_name[index]);
    } else {
      G4cout << " not defined ";
    }
  }
  G4cout << G4endl;
}

const G4String& G4VDecayChannel::GetNoName() const
{
  return noName;
}

#include "Randomize.hh"
G4double G4VDecayChannel::DynamicalMass(G4double massPDG, G4double width, G4double maxDev ) const
{ 
  if (width<=0.0) return massPDG;
  if (maxDev >rangeMass) maxDev = rangeMass;
  if (maxDev <=-1.*rangeMass) return massPDG;  // can not calculate
 
  G4double x = G4UniformRand()*(maxDev+rangeMass) - rangeMass;
  G4double y = G4UniformRand();
  const size_t MAX_LOOP=10000;
  for (size_t loop_counter=0; loop_counter <MAX_LOOP; ++loop_counter){
    if ( y * (width*width*x*x + massPDG*massPDG*width*width) <= massPDG*massPDG*width*width  ) break;
    x = G4UniformRand()*(maxDev+rangeMass) - rangeMass;
    y = G4UniformRand();
  }
  G4double mass = massPDG + x*width;
  return mass;
}
   
G4bool    G4VDecayChannel::IsOKWithParentMass(G4double parentMass)
{
  G4double sumOfDaughterMassMin=0.0;
  CheckAndFillParent();
  CheckAndFillDaughters();
  // skip one body decay
  if (numberOfDaughters==1) return true;
  
  for (G4int index=0; index < numberOfDaughters;  index++) { 
    sumOfDaughterMassMin += 
      G4MT_daughters_mass[index] -rangeMass*G4MT_daughters_width[index];
  }
  return (parentMass >= sumOfDaughterMassMin); 
}

void  G4VDecayChannel::SetBR(G4double value)
{ 
  rbranch = value; 
  if      (rbranch <0.  ) rbranch = 0.0;
  else if (rbranch >1.0 ) rbranch = 1.0;
}

