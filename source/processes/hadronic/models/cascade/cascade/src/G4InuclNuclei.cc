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
// $Id: G4InuclNuclei.cc 71719 2013-06-21 00:01:54Z mkelsey $
//
// 20100301  M. Kelsey -- Add function to create unphysical nuclei for use
//	     as temporary final-state fragments.
// 20100319  M. Kelsey -- Add information message to makeNuclearFragment().
//	     Use new GetBindingEnergy() function instead of bindingEnergy().
// 20100622  M. Kelsey -- Use local "bindingEnergy()" function to call through.
// 20100627  M. Kelsey -- Test for non-physical fragments and abort job.
// 20100630  M. Kelsey -- Use excitation energy in G4Ions
// 20100714  M. Kelsey -- Use G4DynamicParticle::theDynamicalMass to deal with
//	     excitation energy without instantianting "infinite" G4PartDefns.
// 20100719  M. Kelsey -- Change excitation energy without altering momentum
// 20100906  M. Kelsey -- Add fill() functions to rewrite contents
// 20100910  M. Kelsey -- Add clearExitonConfiguration() to fill() functions
// 20100914  M. Kelsey -- Make printout symmetric with G4InuclElemPart,
//		migrate to integer A and Z
// 20100924  M. Kelsey -- Add constructor to copy G4Fragment input, and output
//		functions to create G4Fragment
// 20110214  M. Kelsey -- Replace integer "model" with enum
// 20110308  M. Kelsey -- Follow new G4Fragment interface for hole types
// 20110427  M. Kelsey -- Remove PDG-code warning
// 20110721  M. Kelsey -- Follow base-class ctor change to pass model directly
// 20110829  M. Kelsey -- Add constructor to copy G4V3DNucleus input
// 20110919  M. Kelsey -- Special case:  Allow fill(A=0,Z=0) to make dummy
// 20110922  M. Kelsey -- Add stream argument to printParticle() => print()
// 20121009  M. Kelsey -- Add report of excitons if non-empty
// 20130314  M. Kelsey -- Use G4IonList typedef for fragment map, encapsulate
//		it in a static function with mutexes.
// 20130620  M. Kelsey -- Address Coverity #37503, check self in op=()
// 20140523  M. Kelsey -- Avoid FPE in setExcitationEnergy() for zero Ekin
// 20150608  M. Kelsey -- Label all while loops as terminating.

#include "G4InuclNuclei.hh"
#include "G4AutoLock.hh"
#include "G4Fragment.hh"
#include "G4HadronicException.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4NucleiProperties.hh"
#include "G4Nucleon.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4V3DNucleus.hh"

#include <assert.h>
#include <sstream>
#include <map>

using namespace G4InuclSpecialFunctions;


// Convert contents from (via constructor) and to G4Fragment

G4InuclNuclei::G4InuclNuclei(const G4Fragment& aFragment,
			     G4InuclParticle::Model model)
  : G4InuclParticle() {
  copy(aFragment, model);
}

void G4InuclNuclei::copy(const G4Fragment& aFragment, Model model) {
  fill(aFragment.GetMomentum()/GeV, aFragment.GetA_asInt(),
       aFragment.GetZ_asInt(), aFragment.GetExcitationEnergy(), model);

  // Exciton configuration must be set by hand
  theExitonConfiguration.protonQuasiParticles = aFragment.GetNumberOfCharged();

  theExitonConfiguration.neutronQuasiParticles =
    aFragment.GetNumberOfParticles() - aFragment.GetNumberOfCharged();

  theExitonConfiguration.protonHoles = aFragment.GetNumberOfChargedHoles();

  theExitonConfiguration.neutronHoles =
    aFragment.GetNumberOfHoles() - theExitonConfiguration.protonHoles;
}


// FIXME:  Should we have a local buffer and return by const-reference instead?
G4Fragment G4InuclNuclei::makeG4Fragment() const {
  G4Fragment frag(getA(), getZ(), getMomentum()*GeV);	// From Bertini units

  // Note:  exciton configuration has to be set piece by piece
  frag.SetNumberOfHoles(theExitonConfiguration.protonHoles
			+ theExitonConfiguration.neutronHoles,
			theExitonConfiguration.protonHoles);

  frag.SetNumberOfExcitedParticle(theExitonConfiguration.protonQuasiParticles 
		  + theExitonConfiguration.neutronQuasiParticles,
		  theExitonConfiguration.protonQuasiParticles);

  return frag;
}

G4InuclNuclei::operator G4Fragment() const {
  return makeG4Fragment();
}


// Convert contents from (via constructor) G4V3DNucleus

G4InuclNuclei::G4InuclNuclei(G4V3DNucleus* a3DNucleus,
			     G4InuclParticle::Model model)
  : G4InuclParticle() {
  copy(a3DNucleus, model);
}

void G4InuclNuclei::copy(G4V3DNucleus* a3DNucleus, Model model) {
  if (!a3DNucleus) return;		// Null pointer means no action

  fill(0., a3DNucleus->GetMassNumber(), a3DNucleus->GetCharge(), 0., model);

  // Convert every hit nucleon into an exciton hole
  if (a3DNucleus->StartLoop()) {
    G4Nucleon* nucl = 0;

    /* Loop checking 08.06.2015 MHK */
    while ((nucl = a3DNucleus->GetNextNucleon())) {
      if (nucl->AreYouHit()) {	// Found previously interacted nucleon
	if (nucl->GetParticleType() == G4Proton::Definition())
	  theExitonConfiguration.protonHoles++;

	if (nucl->GetParticleType() == G4Neutron::Definition())
	  theExitonConfiguration.neutronHoles++;
      }
    }
  }
}


// Overwrite data structure (avoids creating/copying temporaries)

void G4InuclNuclei::fill(const G4LorentzVector& mom, G4int a, G4int z,
			 G4double exc, G4InuclParticle::Model model) {
  setDefinition(makeDefinition(a,z));
  setMomentum(mom);
  setExitationEnergy(exc);
  clearExitonConfiguration();
  setModel(model);
}

void G4InuclNuclei::fill(G4double ekin, G4int a, G4int z, G4double exc,
			 G4InuclParticle::Model model) {
  setDefinition(makeDefinition(a,z));
  setKineticEnergy(ekin);
  setExitationEnergy(exc);
  clearExitonConfiguration();
  setModel(model);
}

void G4InuclNuclei::clear() {
  setDefinition(0);
  clearExitonConfiguration();
  setModel(G4InuclParticle::DefaultModel);
}


// Change excitation energy while keeping momentum vector constant

void G4InuclNuclei::setExitationEnergy(G4double e) {
  G4double ekin = getKineticEnergy();		// Current kinetic energy

  G4double emass = getNucleiMass() + e*MeV/GeV;	// From Bertini to G4 units

  // Safety check -- if zero energy, don't do computation
  G4double ekin_new = (ekin == 0.) ? 0.
    : std::sqrt(emass*emass + ekin*(2.*getMass()+ekin)) - emass;

  setMass(emass);	       // Momentum is computed from mass and Ekin
  setKineticEnergy(ekin_new);
}


// Convert nuclear configuration to standard GEANT4 pointer

// WARNING:  Opposite conventions!  G4InuclNuclei uses (A,Z) everywhere, while
//	  G4ParticleTable::GetIon() uses (Z,A)!

G4ParticleDefinition* G4InuclNuclei::makeDefinition(G4int a, G4int z) {
  // SPECIAL CASE:  (0,0) means create dummy without definition
  if (0 == a && 0 == z) return 0;

  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *pd = pTable->GetIonTable()->GetIon(z, a, 0);

  // SPECIAL CASE:  Non-physical nuclear fragment, for final-state return
  if (!pd) pd = makeNuclearFragment(a,z);

  return pd;		// This could return a null pointer if above fails
}


// Shared buffer of nuclear fragments created below, to avoid memory leaks

namespace {
  static std::map<G4int,G4ParticleDefinition*> fragmentList;
  G4Mutex fragListMutex = G4MUTEX_INITIALIZER;
}

// Creates a non-physical pseudo-nucleus, for return as final-state fragment
// from G4IntraNuclearCascader

G4ParticleDefinition* 
G4InuclNuclei::makeNuclearFragment(G4int a, G4int z) {
  if (a<=0 || z<0 || a<z) {
    G4cerr << " >>> G4InuclNuclei::makeNuclearFragment() called with"
	   << " impossible arguments A=" << a << " Z=" << z << G4endl;
    throw G4HadronicException(__FILE__, __LINE__,
			      "G4InuclNuclei impossible A/Z arguments");
  }

  G4int code = G4IonTable::GetNucleusEncoding(z, a);

  // Use local lookup table (see above) to maintain singletons
  // NOTE:  G4ParticleDefinitions don't need to be explicitly deleted
  //        (see comments in G4IonTable.cc::~G4IonTable)

  G4AutoLock fragListLock(&fragListMutex);
  if (fragmentList.find(code) != fragmentList.end()) return fragmentList[code];
  fragListLock.unlock();

  // Name string follows format in G4IonTable.cc::GetIonName(Z,A,E)
  std::stringstream zstr, astr;
  zstr << z;
  astr << a;
  
  G4String name = "Z" + zstr.str() + "A" + astr.str();
  
  G4double mass = getNucleiMass(a,z) *GeV/MeV;	// From Bertini to GEANT4 units
  
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding Excitation-energy
  
  G4Ions* fragPD = new G4Ions(name,       mass, 0., z*eplus,
  			      0,          +1,   0,
			      0,          0,    0,
			      "nucleus",  0,    a, code,
			      true,	  0.,   0,
			      true, "generic",  0,  0.);
  fragPD->SetAntiPDGEncoding(0);

  fragListLock.lock();		    // Protect before saving new fragment
  return (fragmentList[code] = fragPD);     // Store in table for next lookup
}

G4double G4InuclNuclei::getNucleiMass(G4int a, G4int z, G4double exc) {
  // Simple minded mass calculation use constants in CLHEP (all in MeV)
  G4double mass = G4NucleiProperties::GetNuclearMass(a,z) + exc;

  return mass*MeV/GeV;		// Convert from GEANT4 to Bertini units
}

// Assignment operator for use with std::sort()
G4InuclNuclei& G4InuclNuclei::operator=(const G4InuclNuclei& right) {
  if (this != &right) {
    theExitonConfiguration = right.theExitonConfiguration;
    G4InuclParticle::operator=(right);
  }
  return *this;
}

// Dump particle properties for diagnostics

void G4InuclNuclei::print(std::ostream& os) const {
  G4InuclParticle::print(os);
  os << G4endl << " Nucleus: " << getDefinition()->GetParticleName() 
     << " A " << getA() << " Z " << getZ() << " mass " << getMass()
     << " Eex (MeV) " << getExitationEnergy();

  if (!theExitonConfiguration.empty())
    os << G4endl << "         " << theExitonConfiguration;
}
