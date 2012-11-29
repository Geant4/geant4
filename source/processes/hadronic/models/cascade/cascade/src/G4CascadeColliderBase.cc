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
// 20100714  M. Kelsey -- Move functionality from G4VCascadeCollider, and
//		provide conservation-checking here, with wrapper function
//		and control flag.
// 20100721  M. Kelsey -- Use G4CASCADE_CHECK_ECONS to set default control
//		flag for validations.
// 20100923  M. Kelsey -- Migrate to integer A and Z
// 20100925  M. Kelsey -- Add explosion() interfaces for G4Fragment and for
//		(A,Z,E).  Move implementation to latter.  Add Z==0 condition.
// 20110225  M. Kelsey -- Add setVerboseLevel(), calls through to members

#include "G4CascadeColliderBase.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CollisionOutput.hh"
#include "G4Fragment.hh"
#include "G4InteractionCase.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4ios.hh"

using namespace G4InuclSpecialFunctions;


// Constructor and destructor

G4CascadeColliderBase::G4CascadeColliderBase(const char* name, G4int verbose)
  : G4VCascadeCollider(name, verbose),
#ifdef G4CASCADE_CHECK_ECONS    
    doConservationChecks(true),
#else
    doConservationChecks(false),
#endif
    balance(new G4CascadeCheckBalance(0.001, 0.001, name)) {}

G4CascadeColliderBase::~G4CascadeColliderBase() {
  delete balance;
}

void G4CascadeColliderBase::setVerboseLevel(G4int verbose) {
  G4VCascadeCollider::setVerboseLevel(verbose);
  balance->setVerboseLevel(verbose);
}


// Both bullet and target must be hadrons or photons for this to work

G4bool G4CascadeColliderBase::useEPCollider(G4InuclParticle* bullet, 
					    G4InuclParticle* target) const {
  return (dynamic_cast<G4InuclElementaryParticle*>(bullet) &&
	  dynamic_cast<G4InuclElementaryParticle*>(target));
}


// Decide wether nuclear fragment is candidate for G4BigBanger

G4bool G4CascadeColliderBase::explosion(G4InuclNuclei* target) const {
  return target && explosion(target->getA(), target->getZ(), 
			     target->getExitationEnergy());	// in MeV
}

G4bool G4CascadeColliderBase::explosion(G4Fragment* fragment) const {
  return fragment && explosion(fragment->GetA_asInt(), fragment->GetZ_asInt(),
			       fragment->GetExcitationEnergy());     // in MeV
}

G4bool 
G4CascadeColliderBase::explosion(G4int A, G4int Z,
				 G4double excitation) const {
  if (verboseLevel) G4cout << " >>> " << theName << "::explosion ?" << G4endl;

  const G4int a_cut = 20;
  const G4double be_cut = 3.0;

  // Neutron balls, or small fragments with high excitations can explode
  return ((A <= a_cut || Z==0) && 
	  (excitation >= be_cut * bindingEnergy(A,Z))
	  );
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
  if (!doConservationChecks) return true;	// Skip checks if requested

  if (verboseLevel > 1)
    G4cout << " >>> " << theName << "::validateOutput" << G4endl;

  // Show final state particles
  if (verboseLevel > 2) output.printCollisionOutput();

  balance->setVerboseLevel(verboseLevel);
  balance->collide(bullet, target, output);
  return balance->okay();			// Returns false if violations
}

G4bool G4CascadeColliderBase::validateOutput(G4InuclParticle* bullet,
					     G4InuclParticle* target,
		     const std::vector<G4InuclElementaryParticle>& particles) {
  if (!doConservationChecks) return true;	// Skip checks if requested

  if (verboseLevel > 1)
    G4cout << " >>> " << theName << "::validateOutput" << G4endl;

  balance->setVerboseLevel(verboseLevel);
  balance->collide(bullet, target, particles);
  return balance->okay();			// Returns false if violations
}

G4bool G4CascadeColliderBase::validateOutput(G4InuclParticle* bullet,
					     G4InuclParticle* target,
		     const std::vector<G4InuclNuclei>& fragments) {
  if (!doConservationChecks) return true;	// Skip checks if requested

  if (verboseLevel > 1)
    G4cout << " >>> " << theName << "::validateOutput" << G4endl;

  balance->setVerboseLevel(verboseLevel);
  balance->collide(bullet, target, fragments);
  return balance->okay();			// Returns false if violations
}
