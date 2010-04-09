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
// $Id: G4BigBanger.cc,v 1.26 2010-04-09 19:33:11 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100301  M. Kelsey -- In generateBangInSCM(), restore old G4CascMom calcs.
//		for (N-1)th outgoing nucleon.
// 20100319  M. Kelsey -- Use new generateWithRandomAngles for theta,phi stuff
// 20100407  M. Kelsey -- Replace std::vector<> returns with data members.

#include "G4BigBanger.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4ParticleLargerEkin.hh"
#include "G4LorentzConvertor.hh"
#include "G4HadTmpUtil.hh"
#include "G4NucleiProperties.hh"
#include <algorithm>

using namespace G4InuclSpecialFunctions;

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4BigBanger::G4BigBanger()
  : verboseLevel(0) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::G4BigBanger" << G4endl;
  }
}

G4CollisionOutput G4BigBanger::collide(G4InuclParticle* /*bullet*/,
				       G4InuclParticle* target) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::collide" << G4endl;
  }

  // primitive explosion model A -> nucleons to prevent too exotic evaporation

  const G4double small_ekin = 1.0e-6;

  G4CollisionOutput output;
  G4LorentzVector totscm;
  G4LorentzVector totlab;

  if (G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target)) {
    G4double A = nuclei_target->getA();
    G4double Z = nuclei_target->getZ();
    G4LorentzVector PEX = nuclei_target->getMomentum();
    G4double EEXS = nuclei_target->getExitationEnergy();
    G4InuclElementaryParticle dummy(small_ekin, 1);
    G4LorentzConvertor toTheNucleiSystemRestFrame;

    toTheNucleiSystemRestFrame.setBullet(dummy.getMomentum(), dummy.getMass());
    toTheNucleiSystemRestFrame.setTarget(PEX, nuclei_target->getMass());
    toTheNucleiSystemRestFrame.toTheTargetRestFrame();

    //    G4double etot = (EEXS - bindingEnergy(A, Z)) * MeV/GeV;  // To Bertini units
    G4double etot = (EEXS - G4NucleiProperties::GetBindingEnergy(G4lrint(A), G4lrint(Z) ) ) * MeV/GeV;  // To Bertini units
    if (etot < 0.0) {
      G4cout << " Negative energy in BigBanger = " << etot << G4endl;
      etot = 0.0;
    }

    if (verboseLevel > 2) {
      G4cout << " BigBanger: target " << G4endl;
      nuclei_target->printParticle(); 
      G4cout << " BigBanger: a " << A << " z " << Z << " eexs " << EEXS << " etot " <<
	etot << " nm " << nuclei_target->getMass() << G4endl;
    }
  
    generateBangInSCM(etot, A, Z, dummy.getParticleMass(1),
		      dummy.getParticleMass(2));

    if (verboseLevel > 2) {
      G4cout << " particles " << particles.size() << G4endl;
      for(G4int i = 0; i < G4int(particles.size()); i++) 
	particles[i].printParticle();
    }
    if(!particles.empty()) { // convert back to Lab
      particleIterator ipart;

      for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
	if (verboseLevel > 2) {
	  totscm += ipart->getMomentum();
	}
	G4LorentzVector mom = 
	  toTheNucleiSystemRestFrame.backToTheLab(ipart->getMomentum());
	ipart->setMomentum(mom); 

	if (verboseLevel > 2) {
	  mom = ipart->getMomentum();
	  totlab += mom;
	}
      };
      std::sort(particles.begin(), particles.end(), G4ParticleLargerEkin());
      if (verboseLevel > 2) {
	G4cout << " In SCM: total outgoing momentum " << G4endl 
	       << " E " << totscm.e() << " px " << totscm.x()
	       << " py " << totscm.y() << " pz " << totscm.z() << G4endl; 
	G4cout << " In Lab: mom cons " << G4endl 
	       << " E " << PEX.e() + 0.001 * EEXS - totlab.e()
	       << " px " << PEX.x() - totlab.x()
	       << " py " << PEX.y() - totlab.y() 
	       << " pz " << PEX.z() - totlab.z() << G4endl; 
      }
    };	
    output.addOutgoingParticles(particles);
  }
  else {
    G4cout << " BigBanger -> try to bang not nuclei " << G4endl;
  }; 

  return output;
}		     

void G4BigBanger::generateBangInSCM(G4double etot, 
				    G4double a, 
				    G4double z, 
				    G4double mp,
				    G4double mn) {

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
  
  if(ia == 1) {
    // abnormal situation
    G4double m = iz > 0 ? mp : mn;
    G4double pmod = std::sqrt((etot + 2.0 * m) * etot);
    G4LorentzVector mom = generateWithRandomAngles(pmod, m);

    G4int knd = iz > 0 ? 1 : 2;

    particles.push_back(G4InuclElementaryParticle(mom, knd, 8)); // modelId included

    return;
  };  
     
  generateMomentumModules(etot, a, z, mp, mn);
  G4bool bad = true;
  G4int itry = 0;

  while(bad && itry < itry_max) {
    itry++;
    std::vector<G4LorentzVector> scm_momentums;
    G4LorentzVector tot_mom;

    if(ia == 2) {
      // FIXME:  This isn't actually a correct four-vector, wrong mass!
      G4LorentzVector mom = generateWithRandomAngles(momModules[0]);

      tot_mom += mom;		 

      scm_momentums.push_back(mom);

      G4LorentzVector mom1 = -mom;

      scm_momentums.push_back(mom1);  
      bad = false;
    }
    else {
      for(G4int i = 0; i < ia - 2; i++) {
	// FIXME:  This isn't actually a correct four-vector, wrong mass!
	G4LorentzVector mom = generateWithRandomAngles(momModules[i]);

	tot_mom += mom;		 

	scm_momentums.push_back(mom);
      };

      //                handle last two
      G4double tot_mod = tot_mom.rho(); 
      G4double ct = -0.5 * (tot_mod * tot_mod + momModules[ia - 2] * momModules[ia - 2] -
			    momModules[ia - 1] * momModules[ia - 1]) / tot_mod / momModules[ia - 2];

      if (verboseLevel > 2) {
	G4cout << " ct last " << ct << G4endl;
      }
  
      if(std::fabs(ct) < ang_cut) {
	// FIXME:  These are not actually four-vectors, just three-momenta
	G4LorentzVector mom2 = generateWithFixedTheta(ct, momModules[ia - 2]);
	//       rotate to the normal system
	G4LorentzVector apr = tot_mom/tot_mod;
	G4double a_tr = std::sqrt(apr.x()*apr.x() + apr.y()*apr.y());
	G4LorentzVector mom;
	mom.setX(mom2.z()*apr.x() + ( mom2.x()*apr.y() + mom2.y()*apr.z()*apr.x())/a_tr);
	mom.setY(mom2.z()*apr.y() + (-mom2.x()*apr.x() + mom2.y()*apr.z()*apr.y())/a_tr);
	mom.setZ(mom2.z()*apr.z() - mom2.y()*a_tr);

	scm_momentums.push_back(mom);
	//               and the last one
	G4LorentzVector mom1= - mom - tot_mom;
	scm_momentums.push_back(mom1);  
	bad = false;
      };
    };   
    if(!bad) {
      for(G4int i = 0; i < ia; i++) {
	G4int knd = i < iz ? 1 : 2;

	particles.push_back(G4InuclElementaryParticle(scm_momentums[i], knd, 8));
      };
    };
  };  
  if (verboseLevel > 2) {
    if(itry == itry_max) G4cout << " BigBanger -> can not generate bang " << G4endl;
  }

  return;
}
	   
void G4BigBanger::generateMomentumModules(G4double etot, 
					  G4double a, 
					  G4double z, 
					  G4double mp, 
					  G4double mn) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::generateMomentumModules" << G4endl;
  }

  momModules.clear();		// Reset buffer for filling

  G4int ia = G4int(a + 0.1);
  G4int iz = G4int(z + 0.1);

  G4double xtot = 0.0;
  G4double promax = maxProbability(a);
  
  G4int i;
  for(i = 0; i < ia; i++) { 
    G4double x = generateX(ia, a, promax);

    if (verboseLevel > 2) {
      G4cout << " i " << i << " x " << x << G4endl;
    }
    momModules.push_back(x);
    xtot += x;
  };
  for(i = 0; i < ia; i++) {
    G4double m = i < iz ? mp : mn;

    momModules[i] = momModules[i] * etot / xtot;
    momModules[i] = std::sqrt(momModules[i] * (momModules[i] + 2.0 * m));

    if (verboseLevel > 2) {
      G4cout << " i " << i << " pmod " << momModules[i] << G4endl;
    }
  };

  return;  
}

G4double G4BigBanger::xProbability(G4double x, 
				   G4int ia) const {


  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::xProbability" << G4endl;
  }

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

G4double G4BigBanger::generateX(G4int ia, 
				G4double a, 
				G4double promax) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::generateX" << G4endl;
  }

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
