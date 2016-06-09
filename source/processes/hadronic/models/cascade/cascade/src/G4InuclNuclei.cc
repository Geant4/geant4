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
// $Id: G4InuclNuclei.cc,v 1.23 2010-12-15 07:41:11 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
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
// 20110427  M. Kelsey -- Remove PDG-code warning

#include "G4InuclNuclei.hh"
#include "G4Fragment.hh"
#include "G4HadronicException.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include <assert.h>
#include <sstream>
#include <map>

using namespace G4InuclSpecialFunctions;


// Convert contents from (via constructor) and to G4Fragment

G4InuclNuclei::G4InuclNuclei(const G4Fragment& aFragment, G4int model)
  : G4InuclParticle(makeDefinition(aFragment.GetA_asInt(),
				   aFragment.GetZ_asInt()),
		    aFragment.GetMomentum()/GeV) {	// Bertini units
  setExitationEnergy(aFragment.GetExcitationEnergy());
  setModel(model);

  // Exciton configuration must be set by hand
  theExitonConfiguration.protonQuasiParticles = aFragment.GetNumberOfCharged();

  theExitonConfiguration.neutronQuasiParticles =
    aFragment.GetNumberOfCharged() - aFragment.GetNumberOfCharged();

  // Split hole count evenly between protons and neutrons (arbitrary!)
  theExitonConfiguration.protonHoles = aFragment.GetNumberOfHoles()/2;

  theExitonConfiguration.neutronHoles =
    aFragment.GetNumberOfHoles() - theExitonConfiguration.protonHoles;
}

// FIXME:  Should we have a local buffer and return by const-reference instead?
G4Fragment G4InuclNuclei::makeG4Fragment() const {
  G4Fragment frag(getA(), getZ(), getMomentum()*GeV);	// From Bertini units

  // Note:  exciton configuration has to be set piece by piece
  frag.SetNumberOfHoles(theExitonConfiguration.protonHoles
			+ theExitonConfiguration.neutronHoles);

  frag.SetNumberOfParticles(theExitonConfiguration.protonQuasiParticles 
			    + theExitonConfiguration.neutronQuasiParticles);

  frag.SetNumberOfCharged(theExitonConfiguration.protonQuasiParticles);

  return frag;
}

G4InuclNuclei::operator G4Fragment() const {
  return makeG4Fragment();
}


// Overwrite data structure (avoids creating/copying temporaries)

void G4InuclNuclei::fill(const G4LorentzVector& mom, G4int a, G4int z,
			 G4double exc, G4int model) {
  setDefinition(makeDefinition(a,z));
  setMomentum(mom);
  setExitationEnergy(exc);
  clearExitonConfiguration();
  setModel(model);
}

void G4InuclNuclei::fill(G4double ekin, G4int a, G4int z, G4double exc,
			 G4int model) {
  setDefinition(makeDefinition(a,z));
  setKineticEnergy(ekin);
  setExitationEnergy(exc);
  clearExitonConfiguration();
  setModel(model);
}


// Change excitation energy while keeping momentum vector constant

void G4InuclNuclei::setExitationEnergy(G4double e) {
  G4double ekin = getKineticEnergy();		// Current kinetic energy

  G4double emass = getNucleiMass() + e*MeV/GeV;	// From Bertini to G4 units

  // Directly compute new kinetic energy from old
  G4double ekin_new = std::sqrt(emass*emass + ekin*(2.*getMass()+ekin)) - emass;

  setMass(emass);	       // Momentum is computed from mass and Ekin
  setKineticEnergy(ekin_new);
}


// Convert nuclear configuration to standard GEANT4 pointer

// WARNING:  Opposite conventions!  G4InuclNuclei uses (A,Z) everywhere, while
//	  G4ParticleTable::GetIon() uses (Z,A)!

G4ParticleDefinition* G4InuclNuclei::makeDefinition(G4int a, G4int z) {
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *pd = pTable->GetIon(z, a, 0.);

  // SPECIAL CASE:  Non-physical nuclear fragment, for final-state return
  if (!pd) pd = makeNuclearFragment(a,z);

  return pd;		// This could return a null pointer if above fails
}

// Creates a non-standard excited nucleus

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

  // Use local lookup table (see G4IonTable.hh) to maintain singletons
  // NOTE:  G4ParticleDefinitions don't need to be explicitly deleted
  //        (see comments in G4IonTable.cc::~G4IonTable)

  // If correct nucleus already created return it
  static std::map<G4int, G4ParticleDefinition*> fragmentList;
  if (fragmentList.find(code) != fragmentList.end()) return fragmentList[code];

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

  return (fragmentList[code] = fragPD);     // Store in table for next lookup
}

G4double G4InuclNuclei::getNucleiMass(G4int a, G4int z, G4double exc) {
  // Simple minded mass calculation use constants in CLHEP (all in MeV)
  G4double mass = G4NucleiProperties::GetNuclearMass(a,z) + exc;

  return mass*MeV/GeV;		// Convert from GEANT4 to Bertini units
}

// Assignment operator for use with std::sort()
G4InuclNuclei& G4InuclNuclei::operator=(const G4InuclNuclei& right) {
  theExitonConfiguration = right.theExitonConfiguration;
  G4InuclParticle::operator=(right);
  return *this;
}

// Dump particle properties for diagnostics

void G4InuclNuclei::printParticle() const {
  G4InuclParticle::printParticle();
  G4cout << " Nucleus: " << getDefinition()->GetParticleName() 
	 << " A " << getA() << " Z " << getZ() << " mass " << getMass()
	 << " Eex (MeV) " << getExitationEnergy() << G4endl;
}
