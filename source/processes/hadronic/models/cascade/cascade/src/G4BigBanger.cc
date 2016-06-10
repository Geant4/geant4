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
// $Id: G4BigBanger.cc 71942 2013-06-28 19:08:11Z mkelsey $
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
// 20100630  M. Kelsey -- Just do simple boost for target, instead of using
//		G4LorentzConverter with dummy bullet.
// 20100701  M. Kelsey -- Re-throw momentum list, not just angles!
// 20100714  M. Kelsey -- Move conservation checking to base class
// 20100726  M. Kelsey -- Move std::vector<> buffer to .hh file
// 20100923  M. Kelsey -- Migrate to integer A and Z
// 20110214  M. Kelsey -- Follow G4InuclParticle::Model enumerator migration
// 20110806  M. Kelsey -- Pre-allocate buffers to reduce memory churn
// 20110922  M. Kelsey -- Follow G4InuclParticle::print(ostream&) migration
// 20120608  M. Kelsey -- Fix variable-name "shadowing" compiler warnings.
// 20130622  Inherit from G4CascadeDeexciteBase, move to deExcite() interface
//		with G4Fragment
// 20130924  M. Kelsey -- Replace std::pow with G4Pow::powN() for CPU speed
// 20150608  M. Kelsey -- Label all while loops as terminating.

#include <algorithm>

#include "G4BigBanger.hh"
#include "G4SystemOfUnits.hh"
#include "G4CollisionOutput.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4ParticleLargerEkin.hh"
#include "G4Pow.hh"

using namespace G4InuclSpecialFunctions;

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4BigBanger::G4BigBanger() : G4CascadeDeexciteBase("G4BigBanger") {}

void G4BigBanger::deExcite(const G4Fragment& target,
			   G4CollisionOutput& globalOutput) {
  if (verboseLevel) G4cout << " >>> G4BigBanger::deExcite" << G4endl;

  getTargetData(target);
  G4ThreeVector toTheLabFrame = PEX.boostVector();	// From rest to lab

  // This "should" be difference between E-target and sum of m(nucleons)
  G4double etot = (EEXS - bindingEnergy(A,Z)) * MeV/GeV;  // To Bertini units
  if (etot < 0.0) etot = 0.0;
  
  if (verboseLevel > 2) {
    G4cout << " BigBanger: target\n" << target
	   << "\n etot " << etot << G4endl;
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
      G4cout << particles[i] << G4endl;
  }

  if (particles.empty()) {	// No bang!  Don't know why...
    G4cerr << " >>> G4BigBanger unable to process fragment "
	   << target << G4endl;

    // FIXME:  This will violate baryon number, momentum, energy, etc.
    return;
  }

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
    if (verboseLevel > 2) G4cout << *ipart << G4endl;
  }
  
  std::sort(particles.begin(), particles.end(), G4ParticleLargerEkin());

  validateOutput(target, particles);		// Checks <vector> directly
  
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

  globalOutput.addOutgoingParticles(particles);
}		     

void G4BigBanger::generateBangInSCM(G4double etot, G4int a, G4int z) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::generateBangInSCM" << G4endl;
  }

  const G4double ang_cut = 0.9999;
  const G4int itry_max = 1000;
  
  if (verboseLevel > 2) {
    G4cout << " a " << a << " z " << z << G4endl;
  }

  particles.clear();	// Reset output vector before filling
  
  if (a == 1) {		// Special -- bare nucleon doesn't really "explode"
    G4int knd = (z>0) ? 1 : 2;
    particles.push_back(G4InuclElementaryParticle(knd)); // zero momentum
    return;
  }
     
  // NOTE:  If distribution fails, need to regenerate magnitudes and angles!
  //*** generateMomentumModules(etot, a, z);

  scm_momentums.reserve(a);
  G4LorentzVector tot_mom;

  G4bool bad = true;
  G4int itry = 0;
  while(bad && itry < itry_max) {	/* Loop checking 08.06.2015 MHK */
    itry++;
    scm_momentums.clear();

    generateMomentumModules(etot, a, z);
    if (a == 2) {
      // This is only a three-vector, not a four-vector
      G4LorentzVector mom = generateWithRandomAngles(momModules[0]);
      scm_momentums.push_back(mom);
      scm_momentums.push_back(-mom);	// Only safe since three-vector!
      bad = false;
    } else {
      tot_mom *= 0.;		// Easy way to reset accumulator

      for(G4int i = 0; i < a-2; i++) {		// All but last two are thrown
      // This is only a three-vector, not a four-vector
	G4LorentzVector mom = generateWithRandomAngles(momModules[i]);
	scm_momentums.push_back(mom);
	tot_mom += mom;		 
      };

      //                handle last two
      G4double tot_mod = tot_mom.rho(); 
      G4double ct = -0.5*(tot_mod*tot_mod + momModules[a-2]*momModules[a-2]
			  - momModules[a-1]*momModules[a-1]) / tot_mod
	/ momModules[a-2];

      if (verboseLevel > 2) G4cout << " ct last " << ct << G4endl;
  
      if(std::fabs(ct) < ang_cut) {
	// This is only a three-vector, not a four-vector
	G4LorentzVector mom2 = generateWithFixedTheta(ct, momModules[a - 2]);

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
    }	// (a > 2)
  }	// while (bad && itry<itry_max)

  if (!bad) {
    particles.resize(a);	// Use assignment to avoid temporaries
    for(G4int i = 0; i < a; i++) {
      G4int knd = i < z ? 1 : 2;
      particles[i].fill(scm_momentums[i], knd, G4InuclParticle::BigBanger);
    };
  };

  if (verboseLevel > 2) {
    if (itry == itry_max) G4cout << " BigBanger -> can not generate bang " << G4endl;
  }

  return;
}
	   
void G4BigBanger::generateMomentumModules(G4double etot, G4int a, G4int z) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::generateMomentumModules" << G4endl;
  }

  // Proton and neutron masses
  const G4double mp = G4InuclElementaryParticle::getParticleMass(1);
  const G4double mn = G4InuclElementaryParticle::getParticleMass(2);

  momModules.clear();		// Reset buffer for filling

  G4double xtot = 0.0;

  if (a > 2) {			// For "large" nuclei, energy is distributed
    G4double promax = maxProbability(a);

    momModules.resize(a, 0.);	// Pre-allocate to avoid memory churn
    for(G4int i = 0; i < a; i++) { 
      momModules[i] = generateX(a, promax);
      xtot += momModules[i];
      
      if (verboseLevel > 2) {
	G4cout << " i " << i << " x " << momModules[i] << G4endl;
      }
    }
  } else {			// Two-body case is special, must be 50%
    xtot = 1.;
    momModules.push_back(0.5);
    momModules.push_back(0.5);
  }

  for(G4int i = 0; i < a; i++) {
    G4double mass = i < z ? mp : mn;

    momModules[i] *= etot/xtot;
    momModules[i] = std::sqrt(momModules[i] * (momModules[i] + 2.0 * mass));

    if (verboseLevel > 2) {
      G4cout << " i " << i << " pmod " << momModules[i] << G4endl;
    }
  };

  return;  
}

G4double G4BigBanger::xProbability(G4double x, G4int a) const {
  if (verboseLevel > 3) G4cout << " >>> G4BigBanger::xProbability" << G4endl;

  G4Pow* theG4Pow = G4Pow::GetInstance();	// For convenience

  G4double ekpr = 0.0;
  if(x < 1.0 || x > 0.0) {
    ekpr = x * x;

    if (a%2 == 0) { // even A
      ekpr *= std::sqrt(1.0 - x) * theG4Pow->powN((1.0 - x), (3*a-6)/2); 
    }
    else {
      ekpr *= theG4Pow->powN((1.0 - x), (3*a-5)/2);
    };
  }; 
  
  return ekpr;
}

G4double G4BigBanger::maxProbability(G4int a) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::maxProbability" << G4endl;
  }

  return xProbability(2./3./(a-1.0), a);
}

G4double G4BigBanger::generateX(G4int a, G4double promax) const {
  if (verboseLevel > 3) G4cout << " >>> G4BigBanger::generateX" << G4endl;

  const G4int itry_max = 1000;
  G4int itry = 0;
  G4double x;
  
  while(itry < itry_max) {	/* Loop checking 08.06.2015 MHK */
    itry++;
    x = inuclRndm();

    if(xProbability(x, a) >= promax * inuclRndm()) return x;
  };
  if (verboseLevel > 2) {
    G4cout << " BigBanger -> can not generate x " << G4endl;
  }

  return maxProbability(a);
}
