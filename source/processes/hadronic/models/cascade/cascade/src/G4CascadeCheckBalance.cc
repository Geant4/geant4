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
// $Id: G4CascadeCheckBalance.cc,v 1.19 2010-10-19 19:48:32 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Verify and report four-momentum conservation for collision output; uses
// same interface as collision generators.
//
// 20100624  M. Kelsey -- Add baryon number, charge, and kinetic energy
// 20100624  M. Kelsey -- Bug fix:  All checks should be |delta| !
// 20100627  M. Kelsey -- Always report violations on cerr, regardless of
//		verbosity.
// 20100628  M. Kelsey -- Add interface to take list of particles directly,
//		bug fix reporting of charge conservation error.
// 20100630  M. Kelsey -- for nuclei, include excitation energies in total.
// 20100701  M. Kelsey -- Undo previous change, handled by G4InuclNuclei.
// 20100711  M. Kelsey -- Use name of parent collider for reporting messages
// 20100713  M. kelsey -- Hide conservation errors behind verbosity
// 20100715  M. Kelsey -- Use new G4CollisionOutput totals instead of loops,
//		move temporary buffer to be data member
// 20100719  M. Kelsey -- Change zero tolerance to 10 keV instead of 1 keV.
// 20100909  M. Kelsey -- Add interface for both kinds of particle lists
// 20101019  M. Kelsey -- CoVerity report: unitialized constructor

#include "G4CascadeCheckBalance.hh"
#include "globals.hh"
#include "G4CascadParticle.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4CollisionOutput.hh"
#include "G4LorentzVector.hh"
#include <vector>


// Constructor sets acceptance limits

const G4double G4CascadeCheckBalance::tolerance = 1e-6;	// How small is zero?

G4CascadeCheckBalance::G4CascadeCheckBalance(G4double relative,
					     G4double absolute,
					     const char* owner)
  : G4VCascadeCollider(owner), relativeLimit(relative),
    absoluteLimit(absolute), initialBaryon(0), finalBaryon(0),
    initialCharge(0), finalCharge(0) {}


// Pseudo-collision just computes input and output four-vectors

void G4CascadeCheckBalance::collide(G4InuclParticle* bullet,
				    G4InuclParticle* target,
				    G4CollisionOutput& output) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeCheckBalance(" << theName << ")::collide"
	   << G4endl;

  initial *= 0.;	// Fast reset; some colliders only have one pointer
  if (bullet) initial += bullet->getMomentum();
  if (target) initial += target->getMomentum();

  // Baryon number and charge must be computed "by hand"
  initialCharge = 0;
  if (bullet) initialCharge += G4int(bullet->getCharge());
  if (target) initialCharge += G4int(target->getCharge());

  G4InuclElementaryParticle* pbullet =
    dynamic_cast<G4InuclElementaryParticle*>(bullet);
  G4InuclElementaryParticle* ptarget =
    dynamic_cast<G4InuclElementaryParticle*>(target);

  G4InuclNuclei* nbullet = dynamic_cast<G4InuclNuclei*>(bullet);
  G4InuclNuclei* ntarget = dynamic_cast<G4InuclNuclei*>(target);

  initialBaryon =
    ((pbullet ? pbullet->baryon() : nbullet ? nbullet->getA() : 0) +
     (ptarget ? ptarget->baryon() : ntarget ? ntarget->getA() : 0) );

  // Final state totals are computed for us
  final = output.getTotalOutputMomentum();
  finalBaryon = output.getTotalBaryonNumber();
  finalCharge = output.getTotalCharge();

  // Report results
  if (verboseLevel > 2) {
    G4cout << " initial px " << initial.px() << " py " << initial.py()
	   << " pz " << initial.pz() << " E " << initial.e()
	   << " baryon " << initialBaryon << " charge " << initialCharge
	   << "\n   final px " << final.px() << " py " << final.py()
	   << " pz " << final.pz() << " E " << final.e()
	   << " baryon " << finalBaryon << " charge " << finalCharge << G4endl;
  }
}

// Take list of output particles directly (e.g., from G4EPCollider internals)

void G4CascadeCheckBalance::collide(G4InuclParticle* bullet, 
				    G4InuclParticle* target,
    const std::vector<G4InuclElementaryParticle>& particles) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeCheckBalance(" << theName << ")::collide(<vector>)"
	   << G4endl;

  tempOutput.reset();			// Buffer for processing
  tempOutput.addOutgoingParticles(particles);
  collide(bullet, target, tempOutput);
}


// Take list of nuclear fragments directly (e.g., from G4Fissioner internals)

void G4CascadeCheckBalance::collide(G4InuclParticle* bullet, 
				    G4InuclParticle* target,
    const std::vector<G4InuclNuclei>& fragments) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeCheckBalance(" << theName << ")::collide(<vector>)"
	   << G4endl;

  tempOutput.reset();			// Buffer for processing
  tempOutput.addOutgoingNuclei(fragments);
  collide(bullet, target, tempOutput);
}


// Take list of "cparticles" (e.g., from G4NucleiModel internals)

void G4CascadeCheckBalance::collide(G4InuclParticle* bullet,
				    G4InuclParticle* target,
		    const std::vector<G4CascadParticle>& particles) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeCheckBalance(" << theName
	   << ")::collide(<cparticles>)" << G4endl;

  tempOutput.reset();			// Buffer for processing
  tempOutput.addOutgoingParticles(particles);
  collide(bullet, target, tempOutput);
}


// Take lists of both G4InuclEP & G4CP (e.g., from G4IntraNucleiCascader)

void G4CascadeCheckBalance::collide(G4InuclParticle* bullet,
				   G4InuclParticle* target,
	     const std::vector<G4InuclElementaryParticle>& particles,
	     const std::vector<G4CascadParticle>& cparticles) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeCheckBalance(" << theName
	   << ")::collide(<EP>,<CP>)" << G4endl;

  tempOutput.reset();			// Buffer for processing
  tempOutput.addOutgoingParticles(particles);
  tempOutput.addOutgoingParticles(cparticles);
  collide(bullet, target, tempOutput);
}


// Compare relative and absolute violations to limits, and report

G4bool G4CascadeCheckBalance::energyOkay() const {
  G4bool relokay = (std::abs(relativeE()) < relativeLimit);
  G4bool absokay = (std::abs(deltaE()) < absoluteLimit);

  if (verboseLevel && (!relokay || !absokay)) {
    G4cerr << theName << ": Energy conservation: relative " << relativeE()
	   << (relokay ? " conserved" : " VIOLATED")
	   << " absolute " << deltaE()
	   << (absokay ? " conserved" : " VIOLATED") << G4endl;
  } else if (verboseLevel > 2) {
    G4cout << theName << ": Energy conservation: relative " << relativeE()
	   << " conserved absolute " << deltaE() << " conserved" << G4endl;
  }

  return (relokay && absokay);
}

G4bool G4CascadeCheckBalance::ekinOkay() const {
  G4bool relokay = (std::abs(relativeKE()) < relativeLimit);
  G4bool absokay = (std::abs(deltaKE()) < absoluteLimit);

  if (verboseLevel && (!relokay || !absokay)) {
    G4cerr << theName << ": Kinetic energy balance: relative "
	   << relativeKE() << (relokay ? " conserved" : " VIOLATED")
	   << " absolute " << deltaKE()
	   << (absokay ? " conserved" : " VIOLATED") << G4endl;
  } else if (verboseLevel > 2) {
    G4cout << theName << ": Kinetic energy balance: relative "
	   << relativeKE() << " conserved absolute " << deltaKE()
	   << " conserved" << G4endl;
  }

  return (relokay && absokay);
}

G4bool G4CascadeCheckBalance::momentumOkay() const {
  G4bool relokay = (std::abs(relativeP()) < relativeLimit);
  G4bool absokay = (std::abs(deltaP()) < absoluteLimit);

  if (verboseLevel && (!relokay || !absokay)) {
    G4cerr << theName << ": Momentum conservation: relative " << relativeP()
	   << (relokay ? " conserved" : " VIOLATED")
	   << " absolute " << deltaP()
	   << (absokay ? " conserved" : " VIOLATED") << G4endl;
  } else if (verboseLevel > 2) {
    G4cout << theName << ": Momentum conservation: relative " << relativeP()
	   << " conserved absolute " << deltaP() << " conserved" << G4endl;
  }

  return (relokay && absokay);
}

G4bool G4CascadeCheckBalance::baryonOkay() const {
  G4bool bokay = (deltaB() == 0);	// Must be perfect!

  if (verboseLevel && !bokay)
    G4cerr << theName << ": Baryon number VIOLATED " << deltaB() << G4endl;

  return bokay;
}

G4bool G4CascadeCheckBalance::chargeOkay() const {
  G4bool qokay = (deltaQ() == 0);	// Must be perfect!

  if (verboseLevel && !qokay)
    G4cerr << theName << ": Charge conservation VIOLATED " << deltaQ()
	   << G4endl;

  return qokay;
}
