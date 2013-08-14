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
// $Id$
//
// Semi-concrete base class for de-excitation modules, analogous to
// G4CascadeColliderBase.
//
// 20130806  M. Kelsey -- Per A. Dotti, move zero vector to file scope to
//		address thread-collision problem.

#include "G4CascadeDeexciteBase.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CascadeParameters.hh"
#include "G4CollisionOutput.hh"
#include "G4Fragment.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4SystemOfUnits.hh"
#include <vector>

using namespace G4InuclSpecialFunctions;


// Constructor and destructor

G4CascadeDeexciteBase::G4CascadeDeexciteBase(const char* name)
  : G4VCascadeDeexcitation(name), balance(0), A(0), Z(0), EEXS(0) {
  if (G4CascadeParameters::checkConservation())
    balance = new G4CascadeCheckBalance(name);
}

G4CascadeDeexciteBase::~G4CascadeDeexciteBase() {
  delete balance;
}

void G4CascadeDeexciteBase::setVerboseLevel(G4int verbose) {
  G4VCascadeDeexcitation::setVerboseLevel(verbose);
  if (balance) balance->setVerboseLevel(verbose);
}


// Copy pertinent information from G4Fragment for modules

void G4CascadeDeexciteBase::getTargetData(const G4Fragment& target) {
  A = target.GetA_asInt();
  Z = target.GetZ_asInt();
  PEX = target.GetMomentum()/GeV;	// Convert from G4 to Bertini units
  EEXS = target.GetExcitationEnergy();
}


// Create (fill) new G4Fragment with proper momentum/energy handling

namespace {
  static const G4LorentzVector zero(0.,0.,0.,0.);  // File scope avoids churn
}

const G4Fragment& 
G4CascadeDeexciteBase::makeFragment(G4int fragA, G4int fragZ, G4double EX) {
  return makeFragment(zero, fragA, fragZ, EX);
}

const G4Fragment& 
G4CascadeDeexciteBase::makeFragment(G4LorentzVector mom, G4int fragA,
				    G4int fragZ, G4double EX) {
  if (verboseLevel>2) {
    G4cout << " >>> " << theName << "::makeFragment " << mom << " " << fragA
	   << " " << fragZ << " " << EX << G4endl;
  }

  // Adjust four-momentum so that mass is nucleus + excitation
  G4double mass =
    G4InuclNuclei::getNucleiMass(fragA,fragZ) + EX/GeV;
  mom.setVectM(mom.vect(), mass);

  // Overwrite previous fragment contents, zeroing out excitons
  aFragment.SetZandA_asInt(fragZ, fragA);
  aFragment.SetMomentum(mom*GeV);			// Bertini uses GeV!
  aFragment.SetNumberOfHoles(0,0);
  aFragment.SetNumberOfExcitedParticle(0,0);

  return aFragment;
}

// Decide wether nuclear fragment is candidate for G4BigBanger

G4bool G4CascadeDeexciteBase::explosion(const G4Fragment& fragment) const {
  return explosion(fragment.GetA_asInt(), fragment.GetZ_asInt(),
		   fragment.GetExcitationEnergy());     // in MeV
}

G4bool G4CascadeDeexciteBase::explosion(G4int fragA, G4int fragZ,
					 G4double excitation) const {
  if (verboseLevel) G4cout << " >>> " << theName << "::explosion ?" << G4endl;

  const G4int a_cut = 20;
  const G4double be_cut = 3.0;

  // Neutron balls, or small fragments with high excitations can explode
  return ((fragA <= a_cut || fragZ==0) && 
	  (excitation >= be_cut * bindingEnergy(fragA,fragZ))
	  );
}


// Validate output for energy, momentum conservation, etc.

G4bool G4CascadeDeexciteBase::validateOutput(const G4Fragment& target,
					     G4CollisionOutput& output) {
  if (!balance) return true;		// Skip checks unless requested

  if (verboseLevel > 1)
    G4cout << " >>> " << theName << "::validateOutput" << G4endl;

  balance->setVerboseLevel(verboseLevel);
  balance->collide(target, output);
  return balance->okay();			// Returns false if violations
}

G4bool G4CascadeDeexciteBase::validateOutput(const G4Fragment& target,
		     const std::vector<G4InuclElementaryParticle>& particles) {
  if (!balance) return true;		// Skip checks unless requested

  if (verboseLevel > 1)
    G4cout << " >>> " << theName << "::validateOutput" << G4endl;

  balance->setVerboseLevel(verboseLevel);
  balance->collide(target, particles);
  return balance->okay();			// Returns false if violations
}

G4bool G4CascadeDeexciteBase::validateOutput(const G4Fragment& target,
		     const std::vector<G4InuclNuclei>& fragments) {
  if (!balance) return true;		// Skip checks unless requested

  if (verboseLevel > 1)
    G4cout << " >>> " << theName << "::validateOutput" << G4endl;

  balance->setVerboseLevel(verboseLevel);
  balance->collide(target, fragments);
  return balance->okay();			// Returns false if violations
}
