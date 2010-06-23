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
// $Id: G4VCascadeCollider.cc,v 1.3 2010-06-23 19:25:36 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100615  M. Kelsey -- Split constructor to have verbose separately
// 20100623  M. Kelsey -- Use old bindingEnergy() wrapper (now returns
//		G4NucleiProperties::GetBindingEnergy()).

#include "G4VCascadeCollider.hh"
#include "G4InteractionCase.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4ios.hh"

using namespace G4InuclSpecialFunctions;


// Constructors handle verbose announcement

G4VCascadeCollider::G4VCascadeCollider(const char* name)
  : theName(name), verboseLevel(0) {
  if (verboseLevel) G4cout << " >>> " << theName << " ctor " << G4endl;
}


G4VCascadeCollider::G4VCascadeCollider(const char* name, G4int verbose)
  : theName(name), verboseLevel(verbose) {
  if (verboseLevel) G4cout << " >>> " << theName << " ctor " << G4endl;
}


// Both bullet and target must be hadrons or leptons for this to work

G4bool G4VCascadeCollider::useEPCollider(G4InuclParticle* bullet, 
					 G4InuclParticle* target) const {
  return (dynamic_cast<G4InuclElementaryParticle*>(bullet) &&
	  dynamic_cast<G4InuclElementaryParticle*>(target));
}


// Decide wether nuclear fragment is candidate for G4BigBanger

G4bool G4VCascadeCollider::explosion(G4InuclNuclei* target) const {
  if (verboseLevel) {
    G4cout << " >>> " << theName << "::explosion" << G4endl;
  }

  const G4double a_cut = 20.0;
  const G4double be_cut = 3.0;

  G4double a = target->getA();
  G4double z = target->getZ();
  G4double eexs = target->getExitationEnergy();

  // Only small fragments with high excitations can explode
  G4bool explo = ((a <= a_cut) && 
		  (eexs >= be_cut * bindingEnergy(a,z))
		  );

  return explo;
}


// Decide whether bullet-target interaction is candidate for cascade

G4bool 
G4VCascadeCollider::inelasticInteractionPossible(G4InuclParticle* bullet,
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
