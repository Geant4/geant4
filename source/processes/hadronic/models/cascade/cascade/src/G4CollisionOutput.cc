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
// $Id: G4CollisionOutput.cc,v 1.38 2010-12-15 07:41:01 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100309  M. Kelsey -- Introduced bug checking i3 for valid tuning pair
// 20100409  M. Kelsey -- Move non-inlinable code here out of .hh file, replace
//		loop over push_back() with block insert().
// 20100418  M. Kelsey -- Add function to boost output lists to lab frame
// 20100520  M. Kelsey -- Add function to rotate Z axis, from G4Casc.Interface
// 20100620  M. Kelsey -- Add some diagnostics in setOnShell, simplify if's
// 20100630  M. Kelsey -- Use "getExcitationEnergyInGeV()" instead of ".001*"
// 20100701  M. Kelsey -- G4InuclNuclei now includes excitation energy as part
//		of the reported mass and four-vector.
// 20100714  M. Kelsey -- Modify setOnShell() to avoid creating particles
//		with negative kinetic energy.
// 20100715  M. Kelsey -- Add total charge and baryon number functions, and a
//		combined "add()" function to put two of these together
// 20100924  M. Kelsey -- Use "OutgoingNuclei" name consistently, replacing
//		old "TargetFragment".  Add new (reusable) G4Fragment buffer 
//		and access functions for initial post-cascade processing.
//		Move implementation of add() to .cc file.
// 20101019  M. Kelsey -- CoVerity report: unitialized constructor

#include "G4CollisionOutput.hh"
#include "G4CascadParticle.hh"
#include "G4ParticleLargerEkin.hh"
#include "G4LorentzConvertor.hh"
#include "G4LorentzRotation.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include <algorithm>

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;
typedef std::vector<G4InuclNuclei>::iterator nucleiIterator;


G4CollisionOutput::G4CollisionOutput()
  : verboseLevel(0), eex_rest(0), on_shell(false) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CollisionOutput::G4CollisionOutput" << G4endl;
}


G4CollisionOutput& G4CollisionOutput::operator=(const G4CollisionOutput& right)
{
  if (this != &right) {
    verboseLevel = right.verboseLevel;
    outgoingParticles = right.outgoingParticles;
    outgoingNuclei = right.outgoingNuclei; 
    theRecoilFragment = right.theRecoilFragment;
    eex_rest = right.eex_rest;
    on_shell = right.on_shell;
  }
  return *this;
}

void G4CollisionOutput::reset() {
  outgoingNuclei.clear();
  outgoingParticles.clear();

  static const G4Fragment emptyFragment;	// Default ctor is all zeros
  theRecoilFragment = emptyFragment;
}


// Merge two complete objects

void G4CollisionOutput::add(const G4CollisionOutput& right) {
  addOutgoingParticles(right.outgoingParticles);
  addOutgoingNuclei(right.outgoingNuclei);
  theRecoilFragment = right.theRecoilFragment;
}


// Append to lists

void G4CollisionOutput::addOutgoingParticles(const std::vector<G4InuclElementaryParticle>& particles) {
  outgoingParticles.insert(outgoingParticles.end(),
			   particles.begin(), particles.end());
}

void G4CollisionOutput::addOutgoingNuclei(const std::vector<G4InuclNuclei>& nuclea) {
  outgoingNuclei.insert(outgoingNuclei.end(), nuclea.begin(), nuclea.end());
}

// These are primarily for G4IntraNucleiCascader internal checks

void G4CollisionOutput::addOutgoingParticle(const G4CascadParticle& cparticle) {
  addOutgoingParticle(cparticle.getParticle());
}

void G4CollisionOutput::addOutgoingParticles(const std::vector<G4CascadParticle>& cparticles) {
  for (unsigned i=0; i<cparticles.size(); i++)
    addOutgoingParticle(cparticles[i].getParticle());
}

// This comes from PreCompound de-excitation, both particles and nuclei

void G4CollisionOutput::addOutgoingParticles(const G4ReactionProductVector* rproducts) {
  if (!rproducts) return;		// Sanity check, no error if null

  G4ReactionProductVector::const_iterator j;
  for (j=rproducts->begin(); j!=rproducts->end(); ++j) {
    G4ParticleDefinition* pd = (*j)->GetDefinition();

    // FIXME:  This is expensive and unnecessary copying!
    G4DynamicParticle aFragment(pd, (*j)->GetMomentum());
    
    // Nucleons and nuclei are jumbled together in the list
    if (G4InuclElementaryParticle::type(pd)) {
      addOutgoingParticle(G4InuclElementaryParticle(aFragment, 9));
    } else {
      addOutgoingNucleus(G4InuclNuclei(aFragment, 9));
    }
  }
}


G4LorentzVector G4CollisionOutput::getTotalOutputMomentum() const {
  if (verboseLevel > 1)
    G4cout << " >>> G4CollisionOutput::getTotalOutputMomentum" << G4endl;

  G4LorentzVector tot_mom;
  G4int i(0);
  for(i=0; i < G4int(outgoingParticles.size()); i++) {
    tot_mom += outgoingParticles[i].getMomentum();
  }
  for(i=0; i < G4int(outgoingNuclei.size()); i++) {
    tot_mom += outgoingNuclei[i].getMomentum();
  }
  return tot_mom;
}

G4int G4CollisionOutput::getTotalCharge() const {
  if (verboseLevel > 1)
    G4cout << " >>> G4CollisionOutput::getTotalCharge" << G4endl;

  G4int charge = 0;
  G4int i(0);
  for(i=0; i < G4int(outgoingParticles.size()); i++) {
    charge += G4int(outgoingParticles[i].getCharge());
  }
  for(i=0; i < G4int(outgoingNuclei.size()); i++) {
    charge += G4int(outgoingNuclei[i].getCharge());
  }
  return charge;
}

G4int G4CollisionOutput::getTotalBaryonNumber() const {
  if (verboseLevel > 1)
    G4cout << " >>> G4CollisionOutput::getTotalBaryonNumber" << G4endl;

  G4int baryon = 0;
  G4int i(0);
  for(i=0; i < G4int(outgoingParticles.size()); i++) {
    baryon += outgoingParticles[i].baryon();
  }
  for(i=0; i < G4int(outgoingNuclei.size()); i++) {
    baryon += G4int(outgoingNuclei[i].getA());
  }
  return baryon;
}


void G4CollisionOutput::printCollisionOutput() const {
  G4cout << " Output: " << G4endl
	 << " Outgoing Particles: " << outgoingParticles.size() << G4endl;

  G4int i(0);
  for(i=0; i < G4int(outgoingParticles.size()); i++)
    outgoingParticles[i].printParticle(); 

  G4cout << " Outgoing Nuclei: " << outgoingNuclei.size() << G4endl;      
  for(i=0; i < G4int(outgoingNuclei.size()); i++)
    outgoingNuclei[i].printParticle();

  if (theRecoilFragment.GetA() > 0) {
    G4cout << theRecoilFragment << G4endl;
  }
}


// Boost particles and fragment to LAB -- "convertor" must already be configured

void G4CollisionOutput::boostToLabFrame(const G4LorentzConvertor& convertor) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CollisionOutput::boostToLabFrame" << G4endl;

  G4bool withReflection = convertor.reflectionNeeded();

  if (!outgoingParticles.empty()) { 
    particleIterator ipart = outgoingParticles.begin();
    for(; ipart != outgoingParticles.end(); ipart++) {
      G4LorentzVector mom = ipart->getMomentum();
      
      if (withReflection) mom.setZ(-mom.z());
      mom = convertor.rotate(mom);
      mom = convertor.backToTheLab(mom);
      ipart->setMomentum(mom); 
    }

    std::sort(outgoingParticles.begin(), outgoingParticles.end(), G4ParticleLargerEkin());
  }
  
  if (!outgoingNuclei.empty()) { 
    nucleiIterator inuc = outgoingNuclei.begin();
    
    for (; inuc != outgoingNuclei.end(); inuc++) {
      G4LorentzVector mom = inuc->getMomentum(); 
      
      if (withReflection) mom.setZ(-mom.z());
      mom = convertor.rotate(mom);
      mom = convertor.backToTheLab(mom);
      inuc->setMomentum(mom);
    }
  }
}


// Apply LorentzRotation to all particles in event

void G4CollisionOutput::rotateEvent(const G4LorentzRotation& rotate) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CollisionOutput::rotateEvent" << G4endl;

  particleIterator ipart = outgoingParticles.begin();
  for(; ipart != outgoingParticles.end(); ipart++)
    ipart->setMomentum(ipart->getMomentum()*=rotate);

  nucleiIterator inuc = outgoingNuclei.begin();
  for (; inuc != outgoingNuclei.end(); inuc++)
    inuc->setMomentum(inuc->getMomentum()*=rotate);
}


void G4CollisionOutput::trivialise(G4InuclParticle* bullet, 
				   G4InuclParticle* target) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CollisionOutput::trivialize" << G4endl;

  if (G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target)) {
    outgoingNuclei.push_back(*nuclei_target);
  } else {
    G4InuclElementaryParticle* particle =
      dynamic_cast<G4InuclElementaryParticle*>(target);
    outgoingParticles.push_back(*particle);
  }

  if (G4InuclNuclei* nuclei_bullet = dynamic_cast<G4InuclNuclei*>(bullet)) {
    outgoingNuclei.push_back(*nuclei_bullet);
  } else {
    G4InuclElementaryParticle* particle =
      dynamic_cast<G4InuclElementaryParticle*>(bullet);
    outgoingParticles.push_back(*particle);
  }
}


void G4CollisionOutput::setOnShell(G4InuclParticle* bullet, 
				   G4InuclParticle* target) {
  if (verboseLevel > 1)
    G4cout << " >>> G4CollisionOutput::setOnShell" << G4endl;

  const G4double accuracy = 0.00001; // momentum concerves at the level of 10 keV

  on_shell = false;
    
  G4LorentzVector ini_mom = bullet->getMomentum();
  G4LorentzVector momt = target->getMomentum();

  ini_mom += momt;
  
  G4LorentzVector out_mom = getTotalOutputMomentum();
  if(verboseLevel > 2){
    G4cout << " bullet momentum = " << ini_mom.e() <<", "<< ini_mom.x() <<", "<< ini_mom.y()<<", "<< ini_mom.z()<<G4endl;
    G4cout << " target momentum = " << momt.e()<<", "<< momt.x()<<", "<< momt.y()<<", "<< momt.z()<<G4endl;
    G4cout << " Fstate momentum = " << out_mom.e()<<", "<< out_mom.x()<<", "<< out_mom.y()<<", "<< out_mom.z()<<G4endl;
  }

  G4LorentzVector mon_non_cons = ini_mom - out_mom;

  G4double pnc = mon_non_cons.rho();
  G4double enc = mon_non_cons.e();

  setRemainingExitationEnergy();       

  if(verboseLevel > 2){
    printCollisionOutput();
    G4cout << " momentum non conservation: " << G4endl
           << " e " << enc << " p " << pnc << G4endl;
    G4cout << " remaining exitation " << eex_rest << G4endl;
  }

  if(std::fabs(enc) <= accuracy && pnc <= accuracy) {
    on_shell = true;
    return;
  }

  // Adjust "last" particle's four-momentum to balance event
  // ONLY adjust particles with sufficient e or p to remain physical!

  if (verboseLevel > 2) G4cout << " re-balancing four-momenta" << G4endl;

  G4int npart = outgoingParticles.size();
  G4int nnuc = outgoingNuclei.size();
  if (npart > 0) {
    for (G4int ip=npart-1; ip>=0; ip--) {
      if (outgoingParticles[ip].getKineticEnergy()+enc > 0.) {
	G4LorentzVector last_mom = outgoingParticles[ip].getMomentum(); 
	last_mom += mon_non_cons;
	outgoingParticles[ip].setMomentum(last_mom);
	break;
      }
    }
  } else if (nnuc > 0) {
    for (G4int in=nnuc-1; in>=0; in--) {
      if (outgoingNuclei[in].getKineticEnergy()+enc > 0.) {
	G4LorentzVector last_mom = outgoingNuclei[in].getMomentum();
	last_mom += mon_non_cons;
	outgoingNuclei[in].setMomentum(last_mom);
	break;
      }
    }
  }

  out_mom = getTotalOutputMomentum();
  mon_non_cons = ini_mom - out_mom;
  pnc = mon_non_cons.rho();
  enc = mon_non_cons.e();

  if(verboseLevel > 2){
    printCollisionOutput();
    G4cout << " momentum non conservation after (1): " << G4endl 
	   << " e " << enc << " p " << pnc << G4endl;
  }

  // Can energy be balanced just with nuclear excitation?
  G4bool need_hard_tuning = true;
  
  G4double encMeV = mon_non_cons.e() / GeV;	// Excitation below is in MeV
  if (outgoingNuclei.size() > 0) {
    for (G4int i=0; i < G4int(outgoingNuclei.size()); i++) {
      G4double eex = outgoingNuclei[i].getExitationEnergy();
      
      if(eex > 0.0 && eex + encMeV >= 0.0) {
	outgoingNuclei[i].setExitationEnergy(eex+encMeV);
	need_hard_tuning = false;
	break;
      }
    }
    if (need_hard_tuning && encMeV > 0.) {
      outgoingNuclei[0].setExitationEnergy(encMeV);
      need_hard_tuning = false;
    }
  }
  
  if (!need_hard_tuning) {
    on_shell = true;
    return;
  }

  // Momentum (hard) tuning required for energy conservation
  if (verboseLevel > 2)
    G4cout << " trying hard (particle-pair) tuning" << G4endl;

  std::pair<std::pair<G4int, G4int>, G4int> tune_par = selectPairToTune(mon_non_cons.e());
  std::pair<G4int, G4int> tune_particles = tune_par.first;
  G4int mom_ind = tune_par.second;
  
  if(verboseLevel > 2) {
    G4cout << " p1 " << tune_particles.first << " p2 " << tune_particles.second
	   << " ind " << mom_ind << G4endl;
  }

  G4bool tuning_possible = 
    (tune_particles.first >= 0 && tune_particles.second >= 0 &&
     mom_ind >= G4LorentzVector::X);

  if (!tuning_possible) {
    if (verboseLevel > 2) G4cout << " tuning impossible " << G4endl;
    return;
  }
    
  G4LorentzVector mom1 = outgoingParticles[tune_particles.first].getMomentum();
  G4LorentzVector mom2 = outgoingParticles[tune_particles.second].getMomentum();
  G4double newE12 = mom1.e() + mom2.e() + mon_non_cons.e();
  G4double R = 0.5 * (newE12 * newE12 + mom2.e() * mom2.e() - mom1.e() * mom1.e()) / newE12;
  G4double Q = -(mom1[mom_ind] + mom2[mom_ind]) / newE12;
  G4double UDQ = 1.0 / (Q * Q - 1.0);
  G4double W = (R * Q + mom2[mom_ind]) * UDQ;
  G4double V = (mom2.e() * mom2.e() - R * R) * UDQ;
  G4double DET = W * W + V;
  
  if (DET < 0.0) {
    if (verboseLevel > 2) G4cout << " DET < 0 " << G4endl;
    return;
  }

  // Tuning allowed only for non-negative determinant
  G4double x1 = -(W + std::sqrt(DET));
  G4double x2 = -(W - std::sqrt(DET));
  
  // choose the appropriate solution
  G4bool xset = false;
  G4double x = 0.0;
  
  if(mon_non_cons.e() > 0.0) { // x has to be > 0.0
    if(x1 > 0.0) {
      if(R + Q * x1 >= 0.0) {
	x = x1;
	xset = true;
      };
    };  
    if(!xset && x2 > 0.0) {
      if(R + Q * x2 >= 0.0) {
	x = x2;
	xset = true;
      };
    };
  } else { 
    if(x1 < 0.0) {
      if(R + Q * x1 >= 0.) {
	x = x1;
	xset = true;
      };
    };  
    if(!xset && x2 < 0.0) {
      if(R + Q * x2 >= 0.0) {
	x = x2;
	xset = true;
      };
    };
  }	// if(mon_non_cons.e() > 0.0)
  
  if(!xset) {
    if(verboseLevel > 2)
      G4cout << " no appropriate solution found " << G4endl;
    return;
  }	// if(xset)

  // retune momentums
  mom1[mom_ind] += x;
  mom2[mom_ind] -= x;
  outgoingParticles[tune_particles.first ].setMomentum(mom1);
  outgoingParticles[tune_particles.second].setMomentum(mom2);
  out_mom = getTotalOutputMomentum();
  std::sort(outgoingParticles.begin(), outgoingParticles.end(), G4ParticleLargerEkin());
  mon_non_cons = ini_mom - out_mom;
  pnc = mon_non_cons.rho();
  enc = mon_non_cons.e();

  on_shell = (std::fabs(enc) < accuracy || pnc < accuracy);

  if(verboseLevel > 2) {
    G4cout << " momentum non conservation tuning: " << G4endl 
	   << " e " << enc << " p " << pnc 
	   << (on_shell?" success":" FAILURE") << G4endl;
  }
}



void G4CollisionOutput::setRemainingExitationEnergy() { 
  eex_rest = 0.0;
  for(G4int i=0; i < G4int(outgoingNuclei.size()); i++) 
    eex_rest += outgoingNuclei[i].getExitationEnergyInGeV();
}


std::pair<std::pair<G4int, G4int>, G4int> 
G4CollisionOutput::selectPairToTune(G4double de) const {
  if (verboseLevel > 2)
    G4cout << " >>> G4CollisionOutput::selectPairToTune" << G4endl;

  std::pair<G4int, G4int> tup(-1, -1);
  G4int i3 = -1; 
  std::pair<std::pair<G4int, G4int>, G4int> tuner(tup, i3);	// Set invalid

  if (outgoingParticles.size() < 2) return tuner;	// Nothing to do

  G4int ibest1 = -1;
  G4int ibest2 = -1;  
  G4double pbest = 0.0;
  G4double pcut = 0.3 * std::sqrt(1.88 * std::fabs(de));
  G4double p1 = 0.0;
  G4double p2;
  
  for (G4int i = 0; i < G4int(outgoingParticles.size())-1; i++) {
    G4LorentzVector mom1 = outgoingParticles[i].getMomentum();
    
    for (G4int j = i+1; j < G4int(outgoingParticles.size()); j++) {
      G4LorentzVector mom2 = outgoingParticles[j].getMomentum();
      
      for (G4int l = G4LorentzVector::X; l<=G4LorentzVector::Z; l++) {
	if (mom1[l]*mom2[l]<0.0) { 
	  if (std::fabs(mom1[l])>pcut && std::fabs(mom2[l])>pcut) {
	    G4double psum = std::fabs(mom1[l]) + std::fabs(mom2[l]);  
	    
	    if(psum > pbest) {
	      ibest1 = i;
	      ibest2 = j;
	      i3 = l;
	      p1 = mom1[l];
	      p2 = mom2[l];
	      pbest = psum;
	    }	// psum > pbest
	  }	// mom1 and mom2 > pcut
	}	// mom1 ~ -mom2
      }	// for (G4int l ...
    }	// for (G4int j ...
  }	// for (G4int i ...

  if (i3 < 0) return tuner;		
  
  tuner.second = i3;		// Momentum axis for tuning
  
  // NOTE: Sign of de determines order for special case of p1==0.
  if (de > 0.0) {
    tuner.first.first  = (p1>0.) ? ibest1 : ibest2;
    tuner.first.second = (p1>0.) ? ibest2 : ibest1;
  } else {
    tuner.first.first  = (p1<0.) ? ibest2 : ibest1;
    tuner.first.second = (p1<0.) ? ibest1 : ibest2;
  }		 
  
  return tuner;
}
