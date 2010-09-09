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
// $Id: G4CascadeRecoilMaker.cc,v 1.1 2010-09-09 19:11:27 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Collects generated cascade data (using Collider::collide() interface)
// and computes the nuclear recoil kinematics needed to balance the event.
//
// 20100909  M. Kelsey -- Inspired by G4CascadeCheckBalance

#include "G4CascadeRecoilMaker.hh"
#include "globals.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CascadParticle.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4CollisionOutput.hh"
#include "G4LorentzVector.hh"
#include <vector>


// Constructor and destructor

G4CascadeRecoilMaker::G4CascadeRecoilMaker()
  : G4VCascadeCollider("G4CascadeRecoilMaker"),
    recoilA(0.), recoilZ(0.), excitationEnergy(0.) {
  balance = new G4CascadeCheckBalance(0., 0., theName);
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

  balance->setVerboseLevel(verboseLevel);
  balance->collide(bullet, target, particles, cparticles);
  makeRecoilFragment();
}


// Used current event configuration to compute recoil
// NOTE:  CheckBalance uses "final-initial", we want "initial-final"

void G4CascadeRecoilMaker::makeRecoilFragment() {
  recoilZ = -(balance->deltaQ());	// Charge "non-conservation"
  recoilA = -(balance->deltaB());	// Baryon "non-conservation"
  recoilMomentum = -(balance->deltaLV());

  // Nuclear excitation energy 
  excitationEnergy =			// Bertini uses MeV for this, not GeV
    (recoilMomentum.m() - G4InuclNuclei::getNucleiMass(recoilA,recoilZ)) * GeV;

  if (goodRecoil()) {
    theRecoilFragment.fill(recoilMomentum, recoilA, recoilZ, excitationEnergy);
  } else {
    if (verboseLevel > 2 && !wholeEvent())
      G4cout << theName << ": event recoil is not a physical nucleus" << G4endl;
  }
}


// Data quality checks
G4bool G4CascadeRecoilMaker::goodRecoil() const {
  return (recoilA>0 && recoilZ>0 && recoilA >= recoilZ &&
	  recoilMomentum.m() >= G4InuclNuclei::getNucleiMass(recoilA,recoilZ));
}

G4bool G4CascadeRecoilMaker::wholeEvent() const {
  return (recoilA==0 && recoilZ==0 && recoilMomentum.rho() < 1e-6 &&
	  std::abs(recoilMomentum.e()) < 1e-6);
}

