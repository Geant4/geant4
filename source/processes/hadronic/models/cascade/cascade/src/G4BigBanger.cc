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
// $Id: G4BigBanger.cc,v 1.32 2010-06-30 23:07:04 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100301  M. Kelsey -- In generateBangInSCM(), restore old G4CascMom calcs.
//		for (N-1)th outgoing nucleon.
// 20100319  M. Kelsey -- Use new generateWithRandomAngles for theta,phi stuff
// 20100407  M. Kelsey -- Replace std::vector<> returns with data members.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class, clean up code
// 20100628  M. Kelsey -- Use old "bindingEnergy" fn as wrapper, add balance
//		checking after bang.
// 20100630  M. Kelsey -- G4LorentzConverter can't handle isolated target.
//		Just do simple boost.

#include "G4BigBanger.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CollisionOutput.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4ParticleLargerEkin.hh"
#include "G4LorentzConvertor.hh"
#include <algorithm>

using namespace G4InuclSpecialFunctions;

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4BigBanger::G4BigBanger() : G4VCascadeCollider("G4BigBanger") {}

void
G4BigBanger::collide(G4InuclParticle* /*bullet*/, G4InuclParticle* target,
		     G4CollisionOutput& output) {

  if (verboseLevel) G4cout << " >>> G4BigBanger::collide" << G4endl;

  // primitive explosion model A -> nucleons to prevent too exotic evaporation

  G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target);
  if (!nuclei_target) {
    G4cerr << " BigBanger -> try to bang not nuclei " << G4endl;
    return;
  }

  G4double A = nuclei_target->getA();
  G4double Z = nuclei_target->getZ();

  G4LorentzVector PEX = nuclei_target->getMomentum();
  G4double EEXS = nuclei_target->getExitationEnergy();

  // FIXME: Excited "nucleus" hass excitation energy included in mass
  G4double mEXS = nuclei_target->getMass() + EEXS/GeV;
  PEX.setVectM(PEX.vect(), mEXS);

  /*****
  const G4double small_ekin = 1.0e-9;		// To create dummy "bullet"
  G4InuclElementaryParticle dummy(small_ekin, 1);
  *****/
  G4ThreeVector toTheLabFrame = PEX.boostVector();	// From rest to lab

  G4CascadeCheckBalance balance(0.005,0.01);	// Second arg is in GeV
  balance.setVerboseLevel(verboseLevel);

  // This "should" be difference between E-target and sum of m(nucleons)
  G4double etot = (EEXS - bindingEnergy(A,Z)) * MeV/GeV;  // To Bertini units
  if (etot < 0.0) etot = 0.0;
  
  if (verboseLevel > 2) {
    G4cout << " BigBanger: target " << G4endl;
    nuclei_target->printParticle(); 
    G4cout << " etot " << etot << G4endl;
  }

  if (verboseLevel > 3) {
    G4LorentzVector PEXrest = PEX;
    PEXrest.boost(-toTheLabFrame);
    G4cout << " target rest frame: px " << PEXrest.px() << " py "
	   << PEXrest.py() << " pz " << PEXrest.pz() << " E " << PEXrest.e()
	   << G4endl;
  }

  generateBangInSCM(etot, A, Z);
  
  if (verboseLevel > 2) {
    G4cout << " particles " << particles.size() << G4endl;
    for(G4int i = 0; i < G4int(particles.size()); i++) 
      particles[i].printParticle();
  }

  if (particles.empty()) return;	// No bang!  Don't know why...

  // convert back to Lab
  G4LorentzVector totscm;
  G4LorentzVector totlab;

  if (verboseLevel > 2) G4cout << " BigBanger: boosting to lab" << G4endl;

  particleIterator ipart;
  for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
    G4LorentzVector mom = ipart->getMomentum();
    if (verboseLevel > 2) totscm += mom;

    mom.boost(toTheLabFrame);
    if (verboseLevel > 2) totlab += mom;

    ipart->setMomentum(mom); 
    if (verboseLevel > 2) ipart->printParticle();
  }
  
  std::sort(particles.begin(), particles.end(), G4ParticleLargerEkin());
  
  balance.collide(0, target, particles);	// Checks <vector> directly
  balance.okay();				// Report violations
  
  if (verboseLevel > 2) {
    G4cout << " In SCM: total outgoing momentum " << G4endl 
	   << " E " << totscm.e() << " px " << totscm.x()
	   << " py " << totscm.y() << " pz " << totscm.z() << G4endl; 
    G4cout << " In Lab: mom cons " << G4endl 
	   << " E " << PEX.e() - totlab.e()	// PEX now includes EEXS
	   << " px " << PEX.x() - totlab.x()
	   << " py " << PEX.y() - totlab.y() 
	   << " pz " << PEX.z() - totlab.z() << G4endl; 
  }

  output.addOutgoingParticles(particles);
}		     

void G4BigBanger::generateBangInSCM(G4double etot, G4double a, G4double z) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::generateBangInSCM" << G4endl;
  }

  const G4double ang_cut = 0.9999;
  const G4int itry_max = 1000;
  
  G4int ia = int(a + 0.1);
  G4int iz = int(z + 0.1);

  if (verboseLevel > 2) {
    G4cout << " ia " << ia << " iz " << iz << G4endl;
  }

  particles.clear();	// Reset output vector before filling
  
  if (ia == 1) {	// Special -- bare nucleon doesn't really "explode"
    G4int knd = (iz>0) ? 1 : 2;
    particles.push_back(G4InuclElementaryParticle(knd, 8)); // zero momentum
    return;
  }
     
  generateMomentumModules(etot, a, z);

  std::vector<G4LorentzVector> scm_momentums(ia);
  G4LorentzVector tot_mom;

  G4bool bad = true;
  G4int itry = 0;
  while(bad && itry < itry_max) {
    itry++;
    scm_momentums.clear();

    if (ia == 2) {
      // This is only a three-vector, not a four-vector
      G4LorentzVector mom = generateWithRandomAngles(momModules[0]);
      scm_momentums.push_back(mom);
      scm_momentums.push_back(-mom);	// Only safe since three-vector!
      bad = false;
    } else {
      tot_mom *= 0.;		// Easy way to reset accumulator

      for(G4int i = 0; i < ia-2; i++) {		// All but last two are thrown
      // This is only a three-vector, not a four-vector
	G4LorentzVector mom = generateWithRandomAngles(momModules[i]);
	scm_momentums.push_back(mom);
	tot_mom += mom;		 
      };

      //                handle last two
      G4double tot_mod = tot_mom.rho(); 
      G4double ct = -0.5 * (tot_mod*tot_mod + momModules[ia-2]*momModules[ia-2] -
			    momModules[ia-1] * momModules[ia-1]) / tot_mod / momModules[ia-2];

      if (verboseLevel > 2) {
	G4cout << " ct last " << ct << G4endl;
      }
  
      if(std::fabs(ct) < ang_cut) {
	// This is only a three-vector, not a four-vector
	G4LorentzVector mom2 = generateWithFixedTheta(ct, momModules[ia - 2]);

	// rotate to the normal system
	G4LorentzVector apr = tot_mom/tot_mod;
	G4double a_tr = std::sqrt(apr.x()*apr.x() + apr.y()*apr.y());
	G4LorentzVector mom;
	mom.setX(mom2.z()*apr.x() + ( mom2.x()*apr.y() + mom2.y()*apr.z()*apr.x())/a_tr);
	mom.setY(mom2.z()*apr.y() + (-mom2.x()*apr.x() + mom2.y()*apr.z()*apr.y())/a_tr);
	mom.setZ(mom2.z()*apr.z() - mom2.y()*a_tr);

	scm_momentums.push_back(mom);

	// and the last one (again, not actually a four-vector!)
	G4LorentzVector mom1 = -mom - tot_mom;

	scm_momentums.push_back(mom1);  
	bad = false;
      }	// if (abs(ct) < ang_cut)
    }	// (ia > 2)
  }	// while (bad && itry<itry_max)

  if (!bad) {
    for(G4int i = 0; i < ia; i++) {
      G4int knd = i < iz ? 1 : 2;
      particles.push_back(G4InuclElementaryParticle(scm_momentums[i], knd, 8));
    };
  };

  if (verboseLevel > 2) {
    if (itry == itry_max) G4cout << " BigBanger -> can not generate bang " << G4endl;
  }

  return;
}
	   
void G4BigBanger::generateMomentumModules(G4double etot, G4double a, 
					  G4double z) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::generateMomentumModules" << G4endl;
  }

  // Proton and neutron masses
  const G4double mp = G4InuclElementaryParticle::getParticleMass(1);
  const G4double mn = G4InuclElementaryParticle::getParticleMass(2);

  momModules.clear();		// Reset buffer for filling

  G4int ia = G4int(a + 0.1);
  G4int iz = G4int(z + 0.1);

  G4double xtot = 0.0;

  if (ia > 2) {			// For "large" nuclei, energy is distributed
    G4double promax = maxProbability(a);
    
    for(G4int i = 0; i < ia; i++) { 
      G4double x = generateX(ia, a, promax);
      
      if (verboseLevel > 2) {
	G4cout << " i " << i << " x " << x << G4endl;
      }
      momModules.push_back(x);
      xtot += x;
    }
  } else {			// Two-body case is special, must be 50%
    xtot = 1.;
    momModules.push_back(0.5);
    momModules.push_back(0.5);
  }

  for(G4int i = 0; i < ia; i++) {
    G4double m = i < iz ? mp : mn;

    momModules[i] = momModules[i] * etot / xtot;
    momModules[i] = std::sqrt(momModules[i] * (momModules[i] + 2.0 * m));

    if (verboseLevel > 2) {
      G4cout << " i " << i << " pmod " << momModules[i] << G4endl;
    }
  };

  return;  
}

G4double G4BigBanger::xProbability(G4double x, G4int ia) const {
  if (verboseLevel > 3) G4cout << " >>> G4BigBanger::xProbability" << G4endl;

  G4int ihalf = ia / 2;
  G4double ekpr = 0.0;

  if(x < 1.0 || x > 0.0) {
    ekpr = x * x;

    if(2 * ihalf == ia) { // even A
      ekpr *= std::sqrt(1.0 - x) * std::pow((1.0 - x), G4int(G4double(3 * ia - 6) / 2.0)); 
    }
    else {
      ekpr *= std::pow((1.0 - x), G4int(G4double(3 * ia - 5) / 2.0));
    };
  }; 
  
  return ekpr;
}

G4double G4BigBanger::maxProbability(G4double a) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::maxProbability" << G4endl;
  }

  return xProbability(1.0 / (a - 1.0) / 1.5, G4int(a + 0.1));
}

G4double G4BigBanger::generateX(G4int ia, G4double a, 
				G4double promax) const {
  if (verboseLevel > 3) G4cout << " >>> G4BigBanger::generateX" << G4endl;

  const G4int itry_max = 1000;
  G4int itry = 0;
  G4double x;
  
  while(itry < itry_max) {
    itry++;
    x = inuclRndm();

    if(xProbability(x, ia) >= promax * inuclRndm()) return x;
  };
  if (verboseLevel > 2) {
    G4cout << " BigBanger -> can not generate x " << G4endl;
  }

  return maxProbability(a);
}
