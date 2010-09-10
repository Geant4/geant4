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
// $Id: G4CascadeRecoilMaker.cc,v 1.2 2010-09-10 18:03:40 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Collects generated cascade data (using Collider::collide() interface)
// and computes the nuclear recoil kinematics needed to balance the event.
//
// 20100909  M. Kelsey -- Inspired by G4CascadeCheckBalance
// 20100909  M. Kelsey -- Move G4IntraNucleiCascader::goodCase() here, add
//		tolerance for "almost zero" excitation energy

#include "G4CascadeRecoilMaker.hh"
#include "globals.hh"
#include "G4CascadParticle.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CollisionOutput.hh"
#include "G4ExitonConfiguration.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4LorentzVector.hh"
#include <vector>

using namespace G4InuclSpecialFunctions;


// Constructor and destructor

G4CascadeRecoilMaker::G4CascadeRecoilMaker(G4double tolerance)
  : G4VCascadeCollider("G4CascadeRecoilMaker"),
    excTolerance(tolerance), inputEkin(0.),
    recoilA(0.), recoilZ(0.), excitationEnergy(0.) {
  balance = new G4CascadeCheckBalance(tolerance, tolerance, theName);
}
  
G4CascadeRecoilMaker::~G4CascadeRecoilMaker() {
  delete balance;
}


// Standard Collider interface (non-const output "buffer")

void G4CascadeRecoilMaker::collide(G4InuclParticle* bullet,
				    G4InuclParticle* target,
				    G4CollisionOutput& output) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeRecoilMaker::collide" << G4endl;

  // Available energy needed for "goodNucleus()" test at end
  inputEkin = bullet ? bullet->getKineticEnergy() : 0.;

  balance->setVerboseLevel(verboseLevel);
  balance->collide(bullet, target, output);
  makeRecoilFragment();
}

// This is for use with G4IntraNucleiCascader

void G4CascadeRecoilMaker::collide(G4InuclParticle* bullet,
				   G4InuclParticle* target,
	     const std::vector<G4InuclElementaryParticle>& particles,
	     const std::vector<G4CascadParticle>& cparticles) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeRecoilMaker::collide(<EP>,<CP>)" << G4endl;

  // Available energy needed for "goodNucleus()" test at end
  inputEkin = bullet ? bullet->getKineticEnergy() : 0.;

  balance->setVerboseLevel(verboseLevel);
  balance->collide(bullet, target, particles, cparticles);
  makeRecoilFragment();
}


// Used current event configuration to compute recoil
// NOTE:  CheckBalance uses "final-initial", we want "initial-final"

void G4CascadeRecoilMaker::makeRecoilFragment() {
  recoilZ = -(balance->deltaQ());		// Charge "non-conservation"
  recoilA = -(balance->deltaB());		// Baryon "non-conservation"
  recoilMomentum = -(balance->deltaLV());

  // Nuclear excitation energy, which may be "slightly" negative
  if (goodFragment()) {
    excitationEnergy =				// Bertini uses MeV for this
      (recoilMomentum.m() - G4InuclNuclei::getNucleiMass(recoilA,recoilZ)) * GeV;
  } else {
    excitationEnergy = 0.;
  }

  if (std::abs(excitationEnergy) < excTolerance) excitationEnergy = 0.;

  if (verboseLevel > 2) {
    G4cout << "  recoil px " << recoilMomentum.px()
	   << " py " << recoilMomentum.py() << " pz " << recoilMomentum.pz()
	   << " E " << recoilMomentum.e() << " baryon " << recoilA
	   << " charge " << recoilZ
	   << "\n  recoil mass " << recoilMomentum.m()
	   << " 'excitation' energy " << excitationEnergy << G4endl;
  }

  if (goodFragment() && !wholeEvent()) {	// Allow for neg. excitation
    theRecoilFragment.fill(recoilMomentum, recoilA, recoilZ, excitationEnergy);
  } else {
    if (verboseLevel > 2 && !wholeEvent())
      G4cout << theName << ": event recoil is not a physical nucleus" << G4endl;
  }
}


// Modify both local data member and nucleus configuraiton

void G4CascadeRecoilMaker::setRecoilExcitation(G4double Eexc) {
  theRecoilFragment.setExitationEnergy(excitationEnergy = Eexc);
}

void G4CascadeRecoilMaker::
setRecoilExcitonConfig(const G4ExitonConfiguration& excConfig) {
  theRecoilFragment.setExitonConfiguration(excConfig);
}


// Data quality checks

G4bool G4CascadeRecoilMaker::goodFragment() const {
  return (recoilA>0 && recoilZ>=0 && recoilA >= recoilZ);
}

G4bool G4CascadeRecoilMaker::goodRecoil() const {
  return (goodFragment() && excitationEnergy > -excTolerance);
}

G4bool G4CascadeRecoilMaker::wholeEvent() const {
  return (recoilA==0 && recoilZ==0 && recoilMomentum.rho() < 1e-6 &&
	  std::abs(recoilMomentum.e()) < 1e-6);
}


// Determine whether desired nuclear fragment is constructable outcome

G4bool G4CascadeRecoilMaker::goodNucleus() const {
  if (verboseLevel > 2) {
    G4cout << " >>> G4CascadeRecoilMaker::goodNucleus" << G4endl;
  }

  const G4double minExcitation = 0.1*keV;
  const G4double reasonableExcitation = 7.0;	// Multiple of binding energy
  const G4double fractionalExcitation = 0.2;	// Fraction of input to excite

  if (!goodRecoil()) return false;		// Not a sensible nucleus

  if (excitationEnergy < -excTolerance) return false;	// Negative mass-diff

  if (excitationEnergy <= minExcitation) return true;	// Effectively zero

  // Maximum possible excitation energy determined by initial energy
  G4double dm = bindingEnergy(recoilA,recoilZ);
  G4double exc_max0z = fractionalExcitation * inputEkin*GeV;
  G4double exc_dm    = reasonableExcitation * dm;
  G4double exc_max = (exc_max0z > exc_dm) ? exc_max0z : exc_dm;
  
  if (verboseLevel > 3) {
    G4cout << " eexs " << excitationEnergy << " max " << exc_max
	   << " dm " << dm << G4endl;
  }
  
  return (excitationEnergy < exc_max);		// Below maximum possible
}
