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
// $Id: G4CascadeColliderBase.cc 71942 2013-06-28 19:08:11Z mkelsey $
//
// 20100714  M. Kelsey -- Move functionality from G4VCascadeCollider, and
//		provide conservation-checking here, with wrapper function
//		and control flag.
// 20100721  M. Kelsey -- Use G4CASCADE_CHECK_ECONS to set default control
//		flag for validations.
// 20100923  M. Kelsey -- Migrate to integer A and Z
// 20100925  M. Kelsey -- Add explosion() interfaces for G4Fragment and for
//		(A,Z,E).  Move implementation to latter.  Add Z==0 condition.
// 20110225  M. Kelsey -- Add setVerboseLevel(), calls through to members
// 20130621  Move doConservationChecks to G4CascadeParameters, check there
//		before instantiating CheckBalance; change explosion() to
//		use reference, add validateOutput() w/G4Fragment
// 20130622  Move fragment-handling functions to G4CascadeDeexciteBase
// 20140930  Change name from "const char*" to "const G4String"

#include "G4CascadeColliderBase.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CascadeParameters.hh"
#include "G4CollisionOutput.hh"
#include "G4InteractionCase.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclSpecialFunctions.hh"
#include <vector>

using namespace G4InuclSpecialFunctions;


// Constructor and destructor

G4CascadeColliderBase::G4CascadeColliderBase(const G4String& name, G4int verbose)
  : G4VCascadeCollider(name, verbose), balance(0) {
  if (G4CascadeParameters::checkConservation())
    balance = new G4CascadeCheckBalance(name);
}

G4CascadeColliderBase::~G4CascadeColliderBase() {
  delete balance;
}

void G4CascadeColliderBase::setVerboseLevel(G4int verbose) {
  G4VCascadeCollider::setVerboseLevel(verbose);
  if (balance) balance->setVerboseLevel(verbose);
}


// Both bullet and target must be hadrons or photons for this to work

G4bool G4CascadeColliderBase::useEPCollider(G4InuclParticle* bullet, 
					    G4InuclParticle* target) const {
  return (dynamic_cast<G4InuclElementaryParticle*>(bullet) &&
	  dynamic_cast<G4InuclElementaryParticle*>(target));
}


// Decide whether bullet-target interaction is candidate for cascade

G4bool 
G4CascadeColliderBase::inelasticInteractionPossible(G4InuclParticle* bullet,
						    G4InuclParticle* target, 
						    G4double ekin) const {
  if (verboseLevel) {
    G4cout << " >>> " << theName << "::inelasticInteractionPossible" << G4endl;
  }

  // If hadron-hadron collision, defer to ElementaryParticleCollider
  if (useEPCollider(bullet, target)) return true;

  // See which one of the two (or both) is a nucleus, get properties
  // FIXME:  Should set a = baryon() for both, but that's not in base
  G4InuclNuclei* nuclei_bullet = dynamic_cast<G4InuclNuclei*>(bullet);
  G4double ab = nuclei_bullet ? nuclei_bullet->getA() : 1;	// FIXME
  G4double zb = nuclei_bullet ? nuclei_bullet->getZ() : bullet->getCharge();
  
  G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target);
  G4double at = nuclei_target ? nuclei_target->getA() : 1;	// FIXME
  G4double zt = nuclei_target ? nuclei_target->getZ() : target->getCharge();
  
  // VCOL (Coulomb barrier) used for testing if elastic collision necessary
  const G4double coeff = 0.001 * 1.2;

  G4double VCOL = coeff * zt * zb / (G4cbrt(at) + G4cbrt(ab)); 
  
  G4bool possible = true;	// Force inelastic; should be (ekin >= VCOL)

  if (verboseLevel > 3) {
    G4cout << " VCOL: " << VCOL << " ekin: " << ekin << " inelastic possible: "
	   << possible << G4endl;
  }

  return possible;
}


// Validate output for energy, momentum conservation, etc.

G4bool G4CascadeColliderBase::validateOutput(G4InuclParticle* bullet,
					     G4InuclParticle* target,
					     G4CollisionOutput& output) {
  if (!balance) return true;		// Skip checks unless requested

  if (verboseLevel > 1)
    G4cout << " >>> " << theName << "::validateOutput" << G4endl;

  // Show final state particles
  if (verboseLevel > 2) output.printCollisionOutput();

  balance->setVerboseLevel(verboseLevel);
  balance->collide(bullet, target, output);
  return balance->okay();			// Returns false if violations
}

G4bool G4CascadeColliderBase::validateOutput(const G4Fragment& fragment,
					     G4CollisionOutput& output) {
  if (!balance) return true;		// Skip checks unless requested

  if (verboseLevel > 1)
    G4cout << " >>> " << theName << "::validateOutput" << G4endl;

  balance->setVerboseLevel(verboseLevel);
  balance->collide(fragment, output);
  return balance->okay();			// Returns false if violations
}

G4bool G4CascadeColliderBase::validateOutput(G4InuclParticle* bullet,
					     G4InuclParticle* target,
		     const std::vector<G4InuclElementaryParticle>& particles) {
  if (!balance) return true;		// Skip checks unless requested

  if (verboseLevel > 1)
    G4cout << " >>> " << theName << "::validateOutput" << G4endl;

  balance->setVerboseLevel(verboseLevel);
  balance->collide(bullet, target, particles);
  return balance->okay();			// Returns false if violations
}
