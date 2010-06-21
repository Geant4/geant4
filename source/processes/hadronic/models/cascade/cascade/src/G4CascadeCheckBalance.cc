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
// $Id: G4CascadeCheckBalance.cc,v 1.1 2010-06-21 03:40:00 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Verify and report four-momentum conservation for collision output; uses
// same interface as collision generators.

#include "G4CascadeCheckBalance.hh"
#include "globals.hh"
#include "G4InuclParticle.hh"
#include "G4CollisionOutput.hh"
#include "G4LorentzVector.hh"


// Constructor sets acceptance limits

G4CascadeCheckBalance::G4CascadeCheckBalance(G4double relative,
					     G4double absolute)
  : G4VCascadeCollider("G4CascadeCheckBalance"),
    relativeLimit(relative), absoluteLimit(absolute) {}


// Pseudo-collision just computes input and output four-vectors

void G4CascadeCheckBalance::collide(G4InuclParticle* bullet,
				    G4InuclParticle* target,
				    G4CollisionOutput& output) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeCheckBalance::collide" << G4endl;

  initial *= 0.;	// Fast reset
  if (bullet) initial += bullet->getMomentum();		// For evaporators
  if (target) initial += target->getMomentum();

  final = output.getTotalOutputMomentum();

  if (verboseLevel > 2) {
    G4cout << " initial px " << initial.px() << " py " << initial.py()
	   << " pz " << initial.pz() << " E " << initial.e() << G4endl
	   << "   final px " << final.px() << " py " << final.py()
	   << " pz " << final.pz() << " E " << final.e() << G4endl;
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
