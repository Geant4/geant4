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
// $Id: G4InuclNuclei.cc,v 1.7 2010-03-20 22:12:38 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100301  M. Kelsey -- Add function to create unphysical nuclei for use
//	     as temporary final-state fragments.
// 20100319  M. Kelsey -- Add information message to makeNuclearFragment().
//	     Use new GetBindingEnergy() function instead of bindingEnergy().

#include "G4InuclNuclei.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4Ions.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4HadTmpUtil.hh"
#include "G4NucleiProperties.hh"
#include <assert.h>
#include <sstream>
#include <map>

using namespace G4InuclSpecialFunctions;


// Convert nuclear configuration to standard GEANT4 pointer

// WARNING:  Opposite conventions!  G4InuclNuclei uses (A,Z) everywhere, while
//	  G4ParticleTable::GetIon() uses (Z,A)!

G4ParticleDefinition* 
G4InuclNuclei::makeDefinition(G4double a, G4double z, G4double exc) {
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *pd = pTable->GetIon(G4int(z), G4int(a), exc);

  // SPECIAL CASE:  Non-physical nuclear fragment, for final-state return
  if (!pd) pd = makeNuclearFragment(a,z,exc);

  return pd;
}

// Creates a non-physical pseudo-nucleus, for return as final-state fragment
// from G4IntraNuclearCascader

G4ParticleDefinition* 
G4InuclNuclei::makeNuclearFragment(G4double a, G4double z, G4double exc) {
  G4int na=G4int(a), nz=G4int(z), nn=na-nz;	// # nucleon, proton, neutron

  // See G4IonTable.hh::GetNucleusEncoding for explanation
  G4int code = ((100+nz)*1000 + na)*10 + (exc>0.)?1:0;

  // Use local lookup table (see G4IonTable.hh) to maintain singletons
  // NOTE:  G4ParticleDefinitions don't need to be explicitly deleted
  //        (see comments in G4IonTable.cc::~G4IonTable)
  static std::map<G4int, G4ParticleDefinition*> fragmentList;

  if (fragmentList.find(code) != fragmentList.end()) return fragmentList[code];

  // Name string follows format in G4IonTable.cc::GetIonName(Z,A,E)
  std::stringstream zstr, astr, estr;
  zstr << nz;
  astr << na;
  estr << G4int(1000*exc+0.5);	// keV in integer form

  G4String name = "Z" + zstr.str() + "A" + astr.str();
  if (exc>0.) name += "["+estr.str()+"]";

  // Simple minded mass calculation use constants in CLHEP (all in MeV)
  G4double mass = nz*proton_mass_c2 + nn*neutron_mass_c2
    + G4NucleiProperties::GetBindingEnergy(G4lrint(a),G4lrint(z)) + exc;

  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding Excitation-energy

  G4cout << " >>> G4InuclNuclei creating temporary fragment for evaporation "
	 << "with non-standard PDGencoding." << G4endl;

  G4Ions* fragPD = new G4Ions(name,       mass, 0., G4double(nz)*eplus,
			      0,          +1,   0,
			      0,          0,    0,
			      "nucleus",  0,    G4lrint(a), code,
			      true,	  0.,   0,
			      false, "generic", 0,  exc);
  fragPD->SetAntiPDGEncoding(0);

  fragmentList[code] = fragPD;		// Store in table for next lookup
  return fragPD;
}

G4double G4InuclNuclei::getNucleiMass(G4double a, G4double z) {
  G4ParticleDefinition* pd = makeDefinition(a,z);
  return pd ? pd->GetPDGMass()*MeV/GeV : 0.;	// From G4 to Bertini units
}

// Assignment operator for use with std::sort()
G4InuclNuclei& G4InuclNuclei::operator=(const G4InuclNuclei& right) {
  exitationEnergy = right.exitationEnergy;
  theExitonConfiguration = right.theExitonConfiguration;
  G4InuclParticle::operator=(right);
  return *this;
}
