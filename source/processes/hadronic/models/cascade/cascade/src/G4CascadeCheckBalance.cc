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
// $Id: G4CascadeCheckBalance.cc,v 1.2 2010-06-24 20:44:24 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Verify and report four-momentum conservation for collision output; uses
// same interface as collision generators.
//
// 20100624  M. Kelsey -- Add baryon number, charge, and kinetic energy

#include "G4CascadeCheckBalance.hh"
#include "globals.hh"
#include "G4HadTmpUtil.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4CollisionOutput.hh"
#include "G4LorentzVector.hh"


// Constructor sets acceptance limits

G4CascadeCheckBalance::G4CascadeCheckBalance(G4double relative,
					     G4double absolute)
  : G4VCascadeCollider("G4CascadeCheckBalance"),
    relativeLimit(relative), absoluteLimit(absolute),
    initialBaryon(0), finalBaryon(0) {}


// Pseudo-collision just computes input and output four-vectors

void G4CascadeCheckBalance::collide(G4InuclParticle* bullet,
				    G4InuclParticle* target,
				    G4CollisionOutput& output) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeCheckBalance::collide" << G4endl;

  initial *= 0.;	// Fast reset; some colliders only have one pointer
  if (bullet) initial += bullet->getMomentum();
  if (target) initial += target->getMomentum();

  final = output.getTotalOutputMomentum();

  // Baryon number and charge must be computed "by hand"
  G4InuclElementaryParticle* pbullet =
    dynamic_cast<G4InuclElementaryParticle*>(bullet);
  G4InuclElementaryParticle* ptarget =
    dynamic_cast<G4InuclElementaryParticle*>(target);

  G4InuclNuclei* nbullet = dynamic_cast<G4InuclNuclei*>(bullet);
  G4InuclNuclei* ntarget = dynamic_cast<G4InuclNuclei*>(target);

  initialBaryon =
    ((pbullet ? pbullet->baryon() : nbullet ? G4lrint(nbullet->getA()) : 0) +
     (ptarget ? ptarget->baryon() : ntarget ? G4lrint(ntarget->getA()) : 0) );

  initialCharge = bullet->getCharge() + target->getCharge();

  const std::vector<G4InuclNuclei>& nout = output.getNucleiFragments();
  const std::vector<G4InuclElementaryParticle>& pout
    = output.getOutgoingParticles();

  finalBaryon = 0;
  finalCharge = 0;
  for (unsigned ip=0; ip<pout.size(); ip++) {
    finalBaryon += pout[ip].baryon();
    finalCharge += pout[ip].getCharge();
  }

  for (unsigned in=0; in<nout.size(); in++) {
    finalBaryon += G4lrint(nout[in].getA());
    finalCharge += G4lrint(nout[in].getZ());
  }

  // Report results
  if (verboseLevel > 2) {
    G4cout << " initial px " << initial.px() << " py " << initial.py()
	   << " pz " << initial.pz() << " E " << initial.e() << G4endl
	   << " initial baryon number " << initialBaryon << G4endl
	   << "   final px " << final.px() << " py " << final.py()
	   << " pz " << final.pz() << " E " << final.e() << G4endl
	   << "   final baryon number " << finalBaryon << G4endl;
  }
}

// Compare relative and absolute violations to limits, and report

G4bool G4CascadeCheckBalance::energyOkay() const {
  G4bool relokay = (relativeE() < relativeLimit);
  G4bool absokay = (deltaE() < absoluteLimit);

  if (verboseLevel > 2)
    G4cout << " Energy conservation: relative " << relativeE()
	   << (relokay ? " conserved" : " VIOLATED")
	   << " absolute " << deltaE()
	   << (absokay ? " conserved" : " VIOLATED") << G4endl;

  return (relokay && absokay);
}

G4bool ekinOkay() const;
G4bool G4CascadeCheckBalance::ekinOkay() const {
  G4bool relokay = (relativeKE() < relativeLimit);
  G4bool absokay = (deltaKE() < absoluteLimit);

  if (verboseLevel > 2)
    G4cout << " Kinetic energy balance: relative " << relativeKE()
	   << (relokay ? " conserved" : " VIOLATED")
	   << " absolute " << deltaKE()
	   << (absokay ? " conserved" : " VIOLATED") << G4endl;

  return (relokay && absokay);
}

G4bool G4CascadeCheckBalance::momentumOkay() const {
  G4bool relokay = (relativeP() < relativeLimit);
  G4bool absokay = (deltaP() < absoluteLimit);

  if (verboseLevel > 2)
    G4cout << " Momentum conservation: relative " << relativeP()
	   << (relokay ? " conserved" : " VIOLATED")
	   << " absolute " << deltaP()
	   << (absokay ? " conserved" : " VIOLATED") << G4endl;

  return (relokay && absokay);
}

G4bool G4CascadeCheckBalance::baryonOkay() const {
  G4bool bokay = (deltaB() == 0);	// Must be perfect!

  if (!bokay && verboseLevel > 2)
    G4cerr << " Baryon number VIOLATED " << deltaB() << G4endl;

  return bokay;
}

G4bool G4CascadeCheckBalance::chargeOkay() const {
  G4bool qokay = (deltaQ() == 0);	// Must be perfect!

  if (!qokay && verboseLevel > 2)
    G4cerr << " Charge conservation VIOLATED " << deltaB() << G4endl;

  return qokay;
}
