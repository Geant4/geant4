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
// $Id: G4LorentzConvertor.cc,v 1.16 2010-01-12 20:40:34 mkelsey Exp $
//
// 20100108  Michael Kelsey -- Use G4LorentzVector internally

#include "G4LorentzConvertor.hh"
#include "G4HadronicException.hh"

const G4double G4LorentzConvertor::small = 1e-10;

G4LorentzConvertor::G4LorentzConvertor() 
  : verboseLevel(4), degenerated(false) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4LorentzConvertor::G4LorentzConvertor" << G4endl;
  }
}

void G4LorentzConvertor::toTheCenterOfMass() {
   
  if (verboseLevel > 3) {
    G4cout << " >>> G4LorentzConvertor::toTheCenterOfMass" << G4endl;
  }

  G4LorentzVector cm4v = target_mom + bullet_mom;

  ecm_tot = cm4v.m();
  velocity = cm4v.boostVector();
  v2 = velocity.mag2();

  gamma = 1.0 / std::sqrt(std::fabs(1.0 - v2));

  G4LorentzVector targ_cm = target_mom;
  targ_cm.boost(-velocity);
  scm_momentum = -targ_cm.vect();	// Motion of CM relative to target
  pscm = scm_momentum.mag();

  if (verboseLevel > 3)
    G4cout << " i 1 pscm(i) " << scm_momentum.x() << G4endl
	   << " i 2 pscm(i) " << scm_momentum.y() << G4endl
	   << " i 3 pscm(i) " << scm_momentum.z() << G4endl;

  G4double pa = scm_momentum.mag2();
  G4double pb = scm_momentum.dot(velocity);

  ga = v2 - pb * pb / pa;

  degenerated = (ga < small);
  if (degenerated) {
    if (verboseLevel > 3) G4cout << " degenerated case " << G4endl; 
    ga = small;
  } else {
    ga = std::sqrt(ga);
  }

  gb = pb / pscm;
  gbpp = gb / pscm;
  gapp = ga * pscm;

  if (verboseLevel > 3) {
    G4double pv = target_mom.vect().dot(velocity);
    G4cout << " ga " << ga << " v2 " << v2 << " pb " << pb
	   << " pb * pb / pa " << pb * pb / pa << " pv " << pv << G4endl;
  }
}

G4CascadeMomentum G4LorentzConvertor::rotate(const G4LorentzVector& mom) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4LorentzConvertor::rotate(G4CascadeMomentum)" << G4endl
	   << " ga " << ga << " gbpp " << gbpp << " gapp " << gapp << G4endl
	   << " degenerated " << degenerated << G4endl
	   << " before rotation: px " << mom.x() << " py " << mom.y()
	   << " pz " << mom.z() << G4endl;
  }

  G4LorentzVector mom_rot = mom;

  if (!degenerated) {
    G4ThreeVector vscm = velocity - gbpp*scm_momentum;
    G4ThreeVector vxcm = scm_momentum.cross(velocity);

    mom_rot.setVect(mom.x()*vscm/ga + mom.y()*vxcm/gapp +
		    mom.z()*scm_momentum.unit() );
  };

  if (verboseLevel > 3) {
    G4cout << " after rotation: px " << mom_rot.x() << " py " << mom_rot.y() <<
      " pz " << mom_rot.z() << G4endl;
  }

  return mom_rot;
}

G4CascadeMomentum G4LorentzConvertor::rotate(const G4LorentzVector& mom1, 
					     const G4LorentzVector& mom) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4LorentzConvertor::rotate(G4CascadeMomentum,G4CascadeMomentum)" << G4endl
	   << " before rotation: px " << mom.x() << " py " << mom.y()
	   << " pz " << mom.z() << G4endl;
  }

  G4LorentzVector mom_rot = mom;

  G4double pp = mom1.vect().mag2();
  G4double pv = mom1.vect().dot(velocity);

  G4double ga1 = v2 - pv * pv / pp;

  if (verboseLevel > 3) {
    G4cout << " ga1 " << ga1 << " small? " << (ga1 <= small) << G4endl;
  }

  if (ga1 > small) {
    ga1 = std::sqrt(ga1);

    G4double gb1 = pv / pp;

    pp = std::sqrt(pp);

    G4double ga1pp = ga1 * pp;

    if (verboseLevel > 3) {
      G4cout << " gb1 " << gb1 << " ga1pp " << ga1pp << G4endl;
    }

    G4ThreeVector vmom1 = velocity - gb1*mom1;
    G4ThreeVector vxm1  = mom1.vect().cross(velocity);

    mom_rot.setVect(mom.x()*vmom1/ga1 + mom.y()*vxm1/ga1pp +
		    mom.z()*mom1.vect().unit() );
  };

  if (verboseLevel > 3) {
    G4cout << " after rotation: px " << mom_rot.x() << " py " << mom_rot.y() <<
      " pz " << mom_rot.z() << G4endl;
  }

  return mom_rot;
}

void G4LorentzConvertor::toTheTargetRestFrame() {
  if (verboseLevel > 3) {
    G4cout << " >>> G4LorentzConvertor::toTheTargetRestFrame" << G4endl;
  }

  gamma = target_mom.e() / target_mom.m();

  velocity = target_mom.boostVector();
  v2 = velocity.mag2();

  G4LorentzVector bull_cm = bullet_mom;
  bull_cm.boost(-velocity);
  scm_momentum = bull_cm.vect();	// Motion of bullet relative to CM
  plab = pscm = scm_momentum.mag();

  if (verboseLevel > 3)
    G4cout << " rf: i 1 pscm(i) " << scm_momentum.x() << G4endl
	   << " rf: i 2 pscm(i) " << scm_momentum.y() << G4endl
	   << " rf: i 3 pscm(i) " << scm_momentum.z() << G4endl;

  G4double pa = scm_momentum.mag2();
  G4double pb = scm_momentum.dot(velocity);

  ga = v2 - pb * pb / pa;

  degenerated = (ga < small);
  if (degenerated) {
    if (verboseLevel > 3) G4cout << " degenerated case " << G4endl; 
    ga = small;
  } else {
    ga = std::sqrt(ga);
  }

  gb = pb / pscm;
  gbpp = gb / pscm;
  gapp = ga * pscm;   
}

G4CascadeMomentum G4LorentzConvertor::backToTheLab(const G4LorentzVector& mom) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4LorentzConvertor::backToTheLab" << G4endl
	   << " at rest: px " << mom.x() << " py " << mom.y() << " pz " << mom.z()
	   << " e " << mom.e() << G4endl
	   << " v2 " << v2 << G4endl;   
  }

  G4LorentzVector mom1 = mom;
  if (v2 > small) mom1.boost(velocity);

  if (verboseLevel > 3) {
    G4cout << " at lab: px " << mom1.x() << " py " << mom1.y()
	   << " pz " << mom1.z() << G4endl;
  }

  return mom1;
}

G4bool G4LorentzConvertor::reflectionNeeded() const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4LorentzConvertor::reflectionNeeded" << G4endl;
  }

  if (v2 < small && !degenerated) 
    throw G4HadronicException(__FILE__, __LINE__, "G4LorentzConvertor::reflectionNeeded - return value undefined");

  return (v2>=small && (!degenerated || scm_momentum.z() < 0.0));
}







