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
// $Id: G4CascadeRecoilMaker.cc 66241 2012-12-13 18:34:42Z gunter $
//
// Collects generated cascade data (using Collider::collide() interface)
// and computes the nuclear recoil kinematics needed to balance the event.
//
// 20100909  M. Kelsey -- Inspired by G4CascadeCheckBalance
// 20100909  M. Kelsey -- Move G4IntraNucleiCascader::goodCase() here, add
//		tolerance for "almost zero" excitation energy
// 20100910  M. Kelsey -- Drop getRecoilFragment() in favor of user calling
//		makeRecoilFragment() with returned non-const pointer.  Drop
//		handling of excitons.
// 20100921  M. Kelsey -- Return G4InuclNuclei using "makeRecoilNuclei()".
//		Repurpose "makeRecoilFragment()" to return G4Fragment.
// 20100924  M. Kelsey -- Remove unusable G4Fragment::SetExcitationEnergy().
//		Add deltaM to compute mass difference, use excitationEnergy
//		to force G4Fragment four-vector to match.
// 20110214  M. Kelsey -- Follow G4InuclParticle::Model enumerator migration
// 20110303  M. Kelsey -- Add diagnostic messages to goodNucleus().
// 20110308  M. Kelsey -- Follow new G4Fragment interface for hole types
// 20110722  M. Kelsey -- For IntraNucleiCascader, take G4CollOut as argument

#include <vector>

#include "G4CascadeRecoilMaker.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4CascadParticle.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CollisionOutput.hh"
#include "G4Fragment.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4LorentzVector.hh"

using namespace G4InuclSpecialFunctions;


// Constructor and destructor

G4CascadeRecoilMaker::G4CascadeRecoilMaker(G4double tolerance)
  : G4VCascadeCollider("G4CascadeRecoilMaker"),
    excTolerance(tolerance), inputEkin(0.),
    recoilA(0), recoilZ(0), excitationEnergy(0.) {
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
  fillRecoil();
}

// This is for use with G4IntraNucleiCascader

void G4CascadeRecoilMaker::collide(G4InuclParticle* bullet,
				   G4InuclParticle* target,
				   G4CollisionOutput& output,
	     const std::vector<G4CascadParticle>& cparticles) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeRecoilMaker::collide(<EP>,<CP>)" << G4endl;

  // Available energy needed for "goodNucleus()" test at end
  inputEkin = bullet ? bullet->getKineticEnergy() : 0.;

  balance->setVerboseLevel(verboseLevel);
  balance->collide(bullet, target, output, cparticles);
  fillRecoil();
}


// Used current event configuration to construct recoil nucleus
// NOTE:  CheckBalance uses "final-initial", we want "initial-final"

void G4CascadeRecoilMaker::fillRecoil() {
  recoilZ = -(balance->deltaQ());		// Charge "non-conservation"
  recoilA = -(balance->deltaB());		// Baryon "non-conservation"
  recoilMomentum = -(balance->deltaLV());

  theExcitons.clear();		// Discard previous exciton configuraiton

  // Bertini uses MeV for excitation energy
  if (!goodFragment()) excitationEnergy = 0.;
  else excitationEnergy = deltaM() * GeV;

  // Allow for very small negative mass difference, and round to zero
  if (std::abs(excitationEnergy) < excTolerance) excitationEnergy = 0.;

  if (verboseLevel > 2) {
    G4cout << "  recoil px " << recoilMomentum.px()
	   << " py " << recoilMomentum.py() << " pz " << recoilMomentum.pz()
	   << " E " << recoilMomentum.e() << " baryon " << recoilA
	   << " charge " << recoilZ
	   << "\n  recoil mass " << recoilMomentum.m()
	   << " 'excitation' energy " << excitationEnergy << G4endl;
  }
}


// Construct physical nucleus from recoil parameters, if reasonable

G4InuclNuclei* 
G4CascadeRecoilMaker::makeRecoilNuclei(G4InuclParticle::Model model) {
  if (verboseLevel > 1) 
    G4cout << " >>> G4CascadeRecoilMaker::makeRecoilNuclei" << G4endl;

  if (!goodRecoil()) {
    if (verboseLevel > 2 && !wholeEvent())
      G4cout << theName << ": event recoil is not a physical nucleus" << G4endl;

    return 0;		// Null pointer means no fragment
  }

  theRecoilNuclei.fill(recoilMomentum, recoilA, recoilZ,
		       excitationEnergy, model);
  theRecoilNuclei.setExitonConfiguration(theExcitons);

  return &theRecoilNuclei;
}


// Construct pre-compound nuclear fragment from recoil parameters

G4Fragment* G4CascadeRecoilMaker::makeRecoilFragment() {
  if (verboseLevel > 1) 
    G4cout << " >>> G4CascadeRecoilMaker::makeRecoilFragment" << G4endl;

  if (!goodRecoil()) {
    if (verboseLevel > 2 && !wholeEvent())
      G4cout << theName << ": event recoil is not a physical nucleus" << G4endl;

    return 0;		// Null pointer means no fragment
  }

  theRecoilFragment.SetZandA_asInt(recoilZ, recoilA);	// Note convention!

  // User may have overridden excitation energy; force four-momentum to match
  G4double fragMass = 
    G4InuclNuclei::getNucleiMass(recoilA,recoilZ) + excitationEnergy/GeV;

  G4LorentzVector fragMom; fragMom.setVectM(recoilMomentum.vect(), fragMass);
  theRecoilFragment.SetMomentum(fragMom*GeV);	// Bertini uses GeV!

  // Note:  exciton configuration has to be set piece by piece
  //		(arguments are Ntotal,Nproton in both cases)
  G4int nholes = theExcitons.protonHoles+theExcitons.neutronHoles;
  theRecoilFragment.SetNumberOfHoles(nholes, theExcitons.protonHoles);

  G4int nexcit = (theExcitons.protonQuasiParticles
		  + theExcitons.neutronQuasiParticles);
  theRecoilFragment.SetNumberOfExcitedParticle(nexcit,
				    theExcitons.protonQuasiParticles);

  return &theRecoilFragment;
}


// Compute raw mass difference from recoil parameters

G4double G4CascadeRecoilMaker::deltaM() const {
  G4double nucMass = G4InuclNuclei::getNucleiMass(recoilA,recoilZ);
  return (recoilMomentum.m() - nucMass);
}


// Data quality checks

G4bool G4CascadeRecoilMaker::goodFragment() const {
  return (recoilA>0 && recoilZ>=0 && recoilA >= recoilZ);
}

G4bool G4CascadeRecoilMaker::goodRecoil() const {
  return (goodFragment() && excitationEnergy > -excTolerance);
}

G4bool G4CascadeRecoilMaker::wholeEvent() const {
  if (verboseLevel > 2) {
    G4cout << " >>> G4CascadeRecoilMaker::wholeEvent:"
	   << " A " << recoilA << " Z " << recoilZ
	   << " P " << recoilMomentum.rho() << " E " << recoilMomentum.e()
	   << "\n wholeEvent returns "
	   << (recoilA==0 && recoilZ==0 && 
	       recoilMomentum.rho() < excTolerance/GeV &&
	       std::abs(recoilMomentum.e()) < excTolerance/GeV) << G4endl;
  }

  return (recoilA==0 && recoilZ==0 && 
	  recoilMomentum.rho() < excTolerance/GeV &&
	  std::abs(recoilMomentum.e()) < excTolerance/GeV);
}

// Determine whether desired nuclear fragment is constructable outcome

G4bool G4CascadeRecoilMaker::goodNucleus() const {
  if (verboseLevel > 2) {
    G4cout << " >>> G4CascadeRecoilMaker::goodNucleus" << G4endl;
  }

  const G4double minExcitation = 0.1*keV;
  const G4double reasonableExcitation = 7.0;	// Multiple of binding energy
  const G4double fractionalExcitation = 0.2;	// Fraction of input to excite

  if (!goodRecoil()) {
    if (verboseLevel>2) {
      if (!goodFragment()) G4cerr << " goodNucleus: invalid A/Z" << G4endl;
      else if (excitationEnergy < -excTolerance) 
	G4cerr << " goodNucleus: negative excitation" << G4endl;
    }
    return false;				// Not a sensible nucleus
  }

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

  if (verboseLevel > 2 && excitationEnergy >= exc_max)
    G4cerr << " goodNucleus: too much excitation" << G4endl;

  return (excitationEnergy < exc_max);		// Below maximum possible
}
