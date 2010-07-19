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
// $Id: G4InuclNuclei.cc,v 1.14 2010-07-19 22:26:28 mkelsey Exp $
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

#include "G4HadronicException.hh"
#include "G4InuclNuclei.hh"
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


// Change excitation energy while keeping momentum vector constant

void G4InuclNuclei::setExitationEnergy(G4double e) {
  G4double emass = getNucleiMass() + e*MeV/GeV;		// From Bertini to G4 units

  static G4LorentzVector mom;		// Local buffer avoids memory churn
  mom = getMomentum();
  mom.setVectM(mom.vect(), emass);	// Same three-momentum, new mass value

  setMass(emass);
  setMomentum(mom);
}


// Convert nuclear configuration to standard GEANT4 pointer

// WARNING:  Opposite conventions!  G4InuclNuclei uses (A,Z) everywhere, while
//	  G4ParticleTable::GetIon() uses (Z,A)!

G4ParticleDefinition* G4InuclNuclei::makeDefinition(G4double a, G4double z) {
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *pd = pTable->GetIon(G4int(z), G4int(a), 0.);

  // SPECIAL CASE:  Non-physical nuclear fragment, for final-state return
  if (!pd) pd = makeNuclearFragment(a,z);

  return pd;		// This could return a null pointer if above fails
}

// Creates a non-standard excited nucleus

// Creates a non-physical pseudo-nucleus, for return as final-state fragment
// from G4IntraNuclearCascader

G4ParticleDefinition* 
G4InuclNuclei::makeNuclearFragment(G4double a, G4double z) {
  G4int na=G4int(a), nz=G4int(z);	// # nucleons and protons

  if (na<=0 || nz<0 || na<nz) {
    G4cerr << " >>> G4InuclNuclei::makeNuclearFragment() called with"
	   << " impossible arguments A=" << a << " Z=" << z << G4endl;
    throw G4HadronicException(__FILE__, __LINE__,
			      "G4InuclNuclei impossible A/Z arguments");
  }

  G4int code = G4IonTable::GetNucleusEncoding(nz, na);

  // Use local lookup table (see G4IonTable.hh) to maintain singletons
  // NOTE:  G4ParticleDefinitions don't need to be explicitly deleted
  //        (see comments in G4IonTable.cc::~G4IonTable)

  // If correct nucleus already created return it
  static std::map<G4int, G4ParticleDefinition*> fragmentList;
  if (fragmentList.find(code) != fragmentList.end()) return fragmentList[code];

  // Name string follows format in G4IonTable.cc::GetIonName(Z,A,E)
  std::stringstream zstr, astr;
  zstr << nz;
  astr << na;
  
  G4String name = "Z" + zstr.str() + "A" + astr.str();
  
  G4double mass = getNucleiMass(a,z) *GeV/MeV;	// From Bertini to GEANT4 units
  
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding Excitation-energy
  
  G4cout << " >>> G4InuclNuclei creating temporary fragment for evaporation "
	 << "with non-standard PDGencoding." << G4endl;
  
  G4Ions* fragPD = new G4Ions(name,       mass, 0., z*eplus,
  			      0,          +1,   0,
			      0,          0,    0,
			      "nucleus",  0,    na, code,
			      true,	  0.,   0,
			      true, "generic",  0,  0.);
  fragPD->SetAntiPDGEncoding(0);

  return (fragmentList[code] = fragPD);     // Store in table for next lookup
}

G4double G4InuclNuclei::getNucleiMass(G4double a, G4double z, G4double exc) {
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
  G4cout << getDefinition()->GetParticleName() 
	 << " A " << getA() << " Z " << getZ() << " mass " << getMass()
	 << " Eex (MeV) " << getExitationEnergy() << G4endl;
  G4InuclParticle::printParticle();
}
