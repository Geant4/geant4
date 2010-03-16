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
// $Id: G4BigBanger.cc,v 1.21 2010-03-16 22:10:26 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly

#include "G4BigBanger.hh"
#include "G4InuclNuclei.hh"
#include "G4ParticleLargerEkin.hh"
#include "G4LorentzConvertor.hh"
#include <algorithm>

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4BigBanger::G4BigBanger()
  : verboseLevel(1) {
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

    G4double etot = (EEXS - bindingEnergy(A, Z)) * GeV/MeV;

    if (verboseLevel > 2) {
      G4cout << " BigBanger: target " << G4endl;
      nuclei_target->printParticle(); 
      G4cout << " BigBanger: a " << A << " z " << Z << " eexs " << EEXS << " etot " <<
	etot << " nm " << nuclei_target->getMass() << G4endl;
    }
  
    std::vector<G4InuclElementaryParticle> particles = 	    
      generateBangInSCM(etot, A, Z, dummy.getParticleMass(1), dummy.getParticleMass(2));

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

std::vector<G4InuclElementaryParticle>  	    
G4BigBanger::generateBangInSCM(G4double etot, 
			       G4double a, 
			       G4double z, 
			       G4double mp,
			       G4double mn) const {

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
  std::vector<G4InuclElementaryParticle> particles;
  
  if(ia == 1) {
    // abnormal situation
    G4double m = iz > 0 ? mp : mn;
    G4double pmod = std::sqrt((etot + 2.0 * m) * etot);
    G4LorentzVector mom;
    std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
    G4double FI = randomPHI();
    G4double Pt = pmod * COS_SIN.second;

    mom.setX(Pt * std::cos(FI));
    mom.setY(Pt * std::sin(FI));
    mom.setZ(Pt * COS_SIN.first);    

    G4int knd = iz > 0 ? 1 : 2;

    //    particles.push_back(G4InuclElementaryParticle(mom, knd));
    particles.push_back(G4InuclElementaryParticle(mom, knd, 8)); // modelId included

    return particles;
  };  
     
  std::vector<G4double> pmod = generateMomentumModules(etot, a, z, mp, mn);
  G4bool bad = true;
  G4int itry = 0;

  while(bad && itry < itry_max) {
    itry++;
    std::vector<G4LorentzVector> scm_momentums;
    G4LorentzVector tot_mom;

    if(ia == 2) {
      G4LorentzVector mom;
      std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
      double FI = randomPHI();
      double Pt = pmod[0] * COS_SIN.second;

      mom.setX(Pt * std::cos(FI));
      mom.setY(Pt * std::sin(FI));
      mom.setZ(Pt * COS_SIN.first);    

      tot_mom += mom;		 

      scm_momentums.push_back(mom);

      G4LorentzVector mom1 = -mom;

      scm_momentums.push_back(mom1);  
      bad = false;
    }
    else {
      for(G4int i = 0; i < ia - 2; i++) {
	G4LorentzVector mom;
	std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
	G4double FI = randomPHI();
	G4double Pt = pmod[i] * COS_SIN.second;

	mom.setX(Pt * std::cos(FI));
	mom.setY(Pt * std::sin(FI));
	mom.setZ(Pt * COS_SIN.first);    

	tot_mom += mom;		 

	scm_momentums.push_back(mom);
      };

      //                handle last two
      G4double tot_mod = tot_mom.rho(); 
      G4double ct = -0.5 * (tot_mod * tot_mod + pmod[ia - 2] * pmod[ia - 2] -
			    pmod[ia - 1] * pmod[ia - 1]) / tot_mod / pmod[ia - 2];

      if (verboseLevel > 2) {
	G4cout << " ct last " << ct << G4endl;
      }
  
      if(std::fabs(ct) < ang_cut) {
	G4LorentzVector mom2 = generateWithFixedTheta(ct, pmod[ia - 2]);
	//       rotate to the normal system
	G4LorentzVector apr = tot_mom/tot_mod;
	G4LorentzVector mom;
	mom.setVect(mom2.vect().cross(apr));
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

	particles.push_back(G4InuclElementaryParticle(scm_momentums[i], knd));
      };
    };
  };  
  if (verboseLevel > 2) {
    if(itry == itry_max) G4cout << " BigBanger -> can not generate bang " << G4endl;
  }

  return particles;
  
}
	   
std::vector<G4double> G4BigBanger::generateMomentumModules(G4double etot, 
							     G4double a, 
							     G4double z, 
							     G4double mp, 
							     G4double mn) const {


  if (verboseLevel > 3) {
    G4cout << " >>> G4BigBanger::generateMomentumModules" << G4endl;
  }

  G4int ia = int(a + 0.1);
  G4int iz = int(z + 0.1);
  std::vector<G4double> pmod;
  G4double xtot = 0.0;
  G4double promax = maxProbability(a);
  
  G4int i;
  for(i = 0; i < ia; i++) { 
    G4double x = generateX(ia, a, promax);

    if (verboseLevel > 2) {
      G4cout << " i " << i << " x " << x << G4endl;
    }
    pmod.push_back(x);
    xtot += x;
  };
  for(i = 0; i < ia; i++) {
    G4double m = i < iz ? mp : mn;

    pmod[i] = pmod[i] * etot / xtot;
    pmod[i] = std::sqrt(pmod[i] * (pmod[i] + 2.0 * m));

    if (verboseLevel > 2) {
      G4cout << " i " << i << " pmod " << pmod[i] << G4endl;
    }
  };

  return pmod;  
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
