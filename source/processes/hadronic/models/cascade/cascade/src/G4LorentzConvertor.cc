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
// $Id: G4LorentzConvertor.cc,v 1.30 2010-12-15 07:41:13 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100108  Michael Kelsey -- Use G4LorentzVector internally
// 20100112  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100308  M. Kelsey -- Bug fix in toTheTargetRestFrame: scm_momentum should
//		be data member, not local.
// 20100409  M. Kelsey -- Protect std::sqrt(ga) against round-off negatives
// 20100519  M. Kelsey -- Add interfaces to pass G4InuclParticles directly
// 20100617  M. Kelsey -- Add more diagnostic messages with multiple levels
// 20100713  M. Kelsey -- All diagnostic messages should be verbose > 1

#include "G4LorentzConvertor.hh"
#include "G4ThreeVector.hh"
#include "G4HadronicException.hh"
#include "G4InuclParticle.hh"


const G4double G4LorentzConvertor::small = 1.0e-10;

G4LorentzConvertor::G4LorentzConvertor() 
  : verboseLevel(0), velocity(0.), gamma(0.), v2(0.), ecm_tot(0.),
    ga(0.), gb(0.), gbpp(0.), gapp(0.), degenerated(false) {
  if (verboseLevel)
    G4cout << " >>> G4LorentzConvertor::G4LorentzConvertor" << G4endl;
}

G4LorentzConvertor::
G4LorentzConvertor(const G4LorentzVector& bmom, G4double bmass, 
		   const G4LorentzVector& tmom, G4double tmass) 
  : verboseLevel(0), velocity(0.), gamma(0.), v2(0.), ecm_tot(0.),
    ga(0.), gb(0.), gbpp(0.), gapp(0.), degenerated(false) {
  setBullet(bmom, bmass);
  setTarget(tmom, tmass);
}

G4LorentzConvertor::
G4LorentzConvertor(const G4InuclParticle* bullet, 
		   const G4InuclParticle* target) 
  : verboseLevel(0), velocity(0.), gamma(0.), v2(0.), ecm_tot(0.),
    ga(0.), gb(0.), gbpp(0.), gapp(0.), degenerated(false) {
  setBullet(bullet);
  setTarget(target);
}

// Extract four-vectors from input particles
void G4LorentzConvertor::setBullet(const G4InuclParticle* bullet) {
  setBullet(bullet->getMomentum());
}

void G4LorentzConvertor::setTarget(const G4InuclParticle* target) {
  setTarget(target->getMomentum());
}


// Boost bullet and target four-vectors into destired frame

void G4LorentzConvertor::toTheCenterOfMass() {
  if (verboseLevel > 2)
    G4cout << " >>> G4LorentzConvertor::toTheCenterOfMass" << G4endl;

  G4LorentzVector cm4v = target_mom + bullet_mom;
  velocity = cm4v.boostVector();
  if (verboseLevel > 3)
    G4cout << " boost " << velocity.x() << " " << velocity.y() << " "
	   << velocity.z() << G4endl;

  // "SCM" is reverse target momentum in the CM frame
  scm_momentum = target_mom;
  scm_momentum.boost(-velocity);
  scm_momentum.setVect(-scm_momentum.vect());

  if (verboseLevel > 3)
    G4cout << " pscm " << scm_momentum.x() << " " << scm_momentum.y()
	   << " " << scm_momentum.z() << G4endl;

  // Compute kinematic quantities for rotate() functions
  v2 = velocity.mag2();
  ecm_tot = cm4v.m();
  gamma = cm4v.e()/ecm_tot;

  G4double pscm = scm_momentum.rho();
  G4double pa   = pscm*pscm;
  G4double pb   = scm_momentum.vect().dot(velocity);

  G4double gasq = v2-pb*pb/pa;
  ga = (gasq > small*small) ? std::sqrt(gasq) : 0.;

  gb = pb / pscm;
  gbpp = gb / pscm;
  gapp = ga * pscm;

  degenerated = (ga < small);
  if (degenerated && verboseLevel > 2) 
    G4cout << " degenerated case (already in CM frame) " << G4endl; 

  if (verboseLevel > 3) {
    G4double pv = target_mom.vect().dot(velocity);
    G4cout << " ga " << ga << " v2 " << v2 << " pb " << pb
	   << " pb * pb / pa " << pb * pb / pa << " pv " << pv << G4endl;
  }
}

void G4LorentzConvertor::toTheTargetRestFrame() {
  if (verboseLevel > 2)
    G4cout << " >>> G4LorentzConvertor::toTheTargetRestFrame" << G4endl;

  velocity = target_mom.boostVector();
  if (verboseLevel > 3)
    G4cout << " boost " << velocity.x() << " " << velocity.y() << " "
	   << velocity.z() << G4endl;

  // "SCM" is bullet momentum in the target's frame
  scm_momentum = bullet_mom;
  scm_momentum.boost(-velocity);

  if (verboseLevel > 3)
    G4cout << " pseudo-pscm " << scm_momentum.x() << " " << scm_momentum.y()
	   << " " << scm_momentum.z() << G4endl;

  // Compute kinematic quantities for rotate() functions
  v2 = velocity.mag2();
  gamma = target_mom.e() / target_mom.m();

  G4double pscm = scm_momentum.rho();
  G4double pa   = pscm*pscm;
  G4double pb   = velocity.dot(scm_momentum.vect());

  G4double gasq = v2-pb*pb/pa;
  ga = (gasq > small*small) ? std::sqrt(gasq) : 0.;

  gb = pb/pscm;
  gbpp = gb/pscm;
  gapp = ga*pscm;

  degenerated = (ga < small);
  if (degenerated && verboseLevel > 2) 
    G4cout << " degenerated case (already in target frame) " << G4endl; 

  if (verboseLevel > 3) {
    G4cout << " ga " << ga << " v2 " << v2 << " pb " << pb
	   << " pb * pb / pa " << pb * pb / pa << G4endl;
  }
}

G4LorentzVector 
G4LorentzConvertor::backToTheLab(const G4LorentzVector& mom) const {
  if (verboseLevel > 2)
    G4cout << " >>> G4LorentzConvertor::backToTheLab" << G4endl;

  if (verboseLevel > 3)
    G4cout << " at rest: px " << mom.x() << " py " << mom.y() << " pz "
	   << mom.z() << " e " << mom.e() << G4endl
	   << " v2 " << v2 << G4endl;

  G4LorentzVector mom1 = mom;
  if (v2 > small) mom1.boost(velocity);

  if (verboseLevel > 3)
    G4cout << " at lab: px " << mom1.x() << " py " << mom1.y() << " pz "
	   << mom1.z() << G4endl;

  return mom1;
}


// Bullet kinematics in target rest frame (LAB frame, usually)

G4double G4LorentzConvertor::getKinEnergyInTheTRS() const {
  if (verboseLevel > 2)
    G4cout << " >>> G4LorentzConvertor::getKinEnergyInTheTRS" << G4endl;

  G4double pv = bullet_mom.vect().dot(target_mom.vect());
  
   G4double ekin_trf =
    (target_mom.e() * bullet_mom.e() - pv) / target_mom.m()
    - bullet_mom.m();
  
  return ekin_trf; 
}

G4double G4LorentzConvertor::getTRSMomentum() const {
  if (verboseLevel > 2)
    G4cout << " >>> G4LorentzConvertor::getTRSMomentum" << G4endl;

  G4LorentzVector bmom = bullet_mom;
  bmom.boost(-target_mom.boostVector());
  return bmom.rho();
}

G4LorentzVector G4LorentzConvertor::rotate(const G4LorentzVector& mom) const {
  if (verboseLevel > 2)
    G4cout << " >>> G4LorentzConvertor::rotate(G4LorentzVector)" << G4endl;

  if (verboseLevel > 3) {
    G4cout << " ga " << ga << " gbpp " << gbpp << " gapp " << gapp << G4endl
	   << " degenerated " << degenerated << G4endl
	   << " before rotation: px " << mom.x() << " py " << mom.y()
	   << " pz " << mom.z() << G4endl;
  }

  G4LorentzVector mom_rot = mom;
  if (!degenerated) {
    if (verboseLevel > 2)
      G4cout << " rotating to align with reference z axis " << G4endl;

    G4ThreeVector vscm = velocity - gbpp*scm_momentum.vect();
    G4ThreeVector vxcm = scm_momentum.vect().cross(velocity);

    mom_rot.setVect(mom.x()*vscm/ga + mom.y()*vxcm/gapp +
		    mom.z()*scm_momentum.vect().unit() );
  };

  if (verboseLevel > 3) {
    G4cout << " after rotation: px " << mom_rot.x() << " py " << mom_rot.y()
	   << " pz " << mom_rot.z() << G4endl;
  }

  return mom_rot;
}

G4LorentzVector G4LorentzConvertor::rotate(const G4LorentzVector& mom1, 
					   const G4LorentzVector& mom) const {
  if (verboseLevel > 2)
    G4cout << " >>> G4LorentzConvertor::rotate(G4LorentzVector,G4LorentzVector)"
	   << G4endl;

  if (verboseLevel > 3) {
    G4cout << " before rotation: px " << mom.x() << " py " << mom.y()
	   << " pz " << mom.z() << G4endl;
  }

  G4double pp = mom1.vect().mag2();
  G4double pv = mom1.vect().dot(velocity);

  G4double ga1 = v2 - pv * pv / pp;
  if (verboseLevel > 3) {
    G4cout << " ga1 " << ga1 << " small? " << (ga1 <= small) << G4endl;
  }

  G4LorentzVector mom_rot = mom;

  if (ga1 > small) {
    ga1 = std::sqrt(ga1);
    G4double gb1 = pv / pp;

    pp = std::sqrt(pp);
    G4double ga1pp = ga1 * pp;

    if (verboseLevel > 3) {
      G4cout << " gb1 " << gb1 << " ga1pp " << ga1pp << G4endl;
    }

    if (verboseLevel > 2)
      G4cout << " rotating to align with first z axis " << G4endl;

    G4ThreeVector vmom1 = velocity - gb1*mom1;
    G4ThreeVector vxm1  = mom1.vect().cross(velocity);

    mom_rot.setVect(mom.x()*vmom1/ga1 + mom.y()*vxm1/ga1pp +
		    mom.z()*mom1.vect().unit() );
  };

  if (verboseLevel > 3) {
    G4cout << " after rotation: px " << mom_rot.x() << " py " << mom_rot.y()
	   << " pz " << mom_rot.z() << G4endl;
  }

  return mom_rot;
}

G4bool G4LorentzConvertor::reflectionNeeded() const {
  if (verboseLevel > 2)
    G4cout << " >>> G4LorentzConvertor::reflectionNeeded (query)" << G4endl;

  if (verboseLevel > 3) {
    G4cout << " v2 = " << v2 << " SCM z = " << scm_momentum.z()
	   << " degenerated? " << degenerated << G4endl;
  }

  if (v2 < small && !degenerated) 
    throw G4HadronicException(__FILE__, __LINE__, "G4LorentzConvertor::reflectionNeeded - return value undefined");

  if (verboseLevel > 2) {
    G4cout << " reflection across XY is"
	   << ((v2>=small && (!degenerated || scm_momentum.z()<0.0))?"":" NOT")
	   << " needed" << G4endl;
  }

  return (v2>=small && (!degenerated || scm_momentum.z()<0.0));
}


// Reporting functions for diagnostics

void G4LorentzConvertor::printBullet() const {
  G4cout << " G4LC bullet: px " << bullet_mom.px() << " py " << bullet_mom.py()
	 << " pz " << bullet_mom.pz() << " e " << bullet_mom.e()
	 << " mass " << bullet_mom.m() << G4endl;
  }

void G4LorentzConvertor::printTarget() const {
  G4cout << " G4LC target: px " << target_mom.px() << " py " << target_mom.py()
	 << " pz " << target_mom.pz() << " e " << target_mom.e()
	 << " mass " << target_mom.m() << G4endl;
}

