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
// $Id: G4CollisionOutput.cc,v 1.19 2010-03-16 22:10:26 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly

#include "G4CollisionOutput.hh"
#include "G4ParticleLargerEkin.hh"
#include <algorithm>

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;


G4CollisionOutput::G4CollisionOutput()
  : verboseLevel(0) {
  
  if (verboseLevel > 3)
    G4cout << " >>> G4CollisionOutput::G4CollisionOutput" << G4endl;
}


G4CollisionOutput& G4CollisionOutput::operator=(const G4CollisionOutput& right)
{
  if (this != &right) {
    verboseLevel = right.verboseLevel;
    outgoingParticles = right.outgoingParticles;
    nucleiFragments = right.nucleiFragments; 
    eex_rest = right.eex_rest;
    on_shell = right.on_shell;
  }
  return *this;
}


G4LorentzVector G4CollisionOutput::getTotalOutputMomentum() const {
  G4LorentzVector tot_mom;
  double eex_r = 0.0;
  G4int i(0);
  for(i = 0; i < G4int(outgoingParticles.size()); i++) {
    tot_mom += outgoingParticles[i].getMomentum();
  }
  for(i = 0; i < G4int(nucleiFragments.size()); i++) {
    tot_mom += nucleiFragments[i].getMomentum();
    eex_r += 0.001 * nucleiFragments[i].getExitationEnergy();
  }
  tot_mom.setE(tot_mom.e() + eex_r);
  return tot_mom;
}

void G4CollisionOutput::setOnShell(G4InuclParticle* bullet, 
				   G4InuclParticle* target) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4CollisionOutput::setOnShell" << G4endl;
  }

  const G4double accuracy = 0.00001; // momentum concerves at the level of 1 eV

  on_shell = false;
    
  G4LorentzVector ini_mom = bullet->getMomentum();
  G4LorentzVector momt = target->getMomentum();

  ini_mom += momt;
  
  G4LorentzVector out_mom = getTotalOutputMomentum();
  G4LorentzVector mon_non_cons;
  if(verboseLevel > 2){
    G4cout << " bullet momentum = " << ini_mom.e() <<", "<< ini_mom.x() <<", "<< ini_mom.y()<<", "<< ini_mom.z()<<G4endl;
    G4cout << " target momentum = " << momt.e()<<", "<< momt.x()<<", "<< momt.y()<<", "<< momt.z()<<G4endl;
    G4cout << " Fstate momentum = " << out_mom.e()<<", "<< out_mom.x()<<", "<< out_mom.y()<<", "<< out_mom.z()<<G4endl;
  }

  mon_non_cons = ini_mom - out_mom;

  G4double pnc = mon_non_cons.rho();

  setRemainingExitationEnergy();       

  if(verboseLevel > 2){
    printCollisionOutput();
    G4cout << " momentum non conservation: " << G4endl
           << " e " << mon_non_cons.e()
	   << " p " << pnc << G4endl;
    G4cout << " remaining exitation " << eex_rest << G4endl;
  }

  if(std::fabs(mon_non_cons.e()) > accuracy || pnc > accuracy) { // renormalization
    G4int npart = outgoingParticles.size();
 
    if(npart > 0) {

      G4LorentzVector last_mom = outgoingParticles[npart - 1].getMomentum(); 

      last_mom += mon_non_cons;
      outgoingParticles[npart - 1].setMomentum(last_mom);
    }
    else {

      G4int nnuc = nucleiFragments.size();
 
      if(nnuc > 0) {

        G4LorentzVector last_mom = nucleiFragments[nnuc - 1].getMomentum();

        last_mom += mon_non_cons;
	nucleiFragments[nnuc - 1].setMomentum(last_mom);
	nucleiFragments[nnuc - 1].setEnergy();
      };
    }; 
    out_mom = getTotalOutputMomentum();
    mon_non_cons = ini_mom - out_mom;
    pnc = mon_non_cons.rho();

    if(verboseLevel > 2){
      printCollisionOutput();
      G4cout << " momentum non conservation after (1): " << G4endl 
	     << " e " << mon_non_cons.e() << " p " << pnc << G4endl;
    }
    G4bool need_hard_tuning = true;
    
    if(nucleiFragments.size() > 0) {
      for(G4int i = 0; i < G4int(nucleiFragments.size()); i++) {

        G4double eex = nucleiFragments[i].getExitationEnergyInGeV();

	if(eex > 0.0 && eex + mon_non_cons.e() >= 0.0) {
	  nucleiFragments[i].setExitationEnergy(1000.0 * (eex + mon_non_cons.e()));
	  need_hard_tuning = false;
	  break;
	};
      };
      if(need_hard_tuning && mon_non_cons.e() > 0.) {
	nucleiFragments[0].setExitationEnergy(1000.0 * mon_non_cons.e());
	need_hard_tuning = false;
      };
    };
    
    if(need_hard_tuning) {
     
      std::pair<std::pair<G4int, G4int>, G4int> tune_par = selectPairToTune(mon_non_cons.e());
      std::pair<G4int, G4int> tune_particles = tune_par.first;
      G4int mom_ind = tune_par.second;

      if(verboseLevel > 2) {
	G4cout << " p1 " << tune_particles.first << " p2 " << tune_particles.second
	       << " ind " << mom_ind << G4endl;
      }
      if(tune_particles.first >= 0 && tune_particles.second >= 0 &&
	 mom_ind >= 1) { // tunning possible

        G4LorentzVector mom1 = outgoingParticles[tune_particles.first].getMomentum();
        G4LorentzVector mom2 = outgoingParticles[tune_particles.second].getMomentum();
        G4double newE12 = mom1.e() + mom2.e() + mon_non_cons.e();
        G4double R = 0.5 * (newE12 * newE12 + mom2.e() * mom2.e() - mom1.e() * mom1.e()) / newE12;
        G4double Q = -(mom1[mom_ind] + mom2[mom_ind]) / newE12;
        G4double UDQ = 1.0 / (Q * Q - 1.0);
        G4double W = (R * Q + mom2[mom_ind]) * UDQ;
        G4double V = (mom2.e() * mom2.e() - R * R) * UDQ;
        G4double DET = W * W + V;

        if(DET > 0.0) {

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
	  }
	  else { 
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
	  };  
          if(xset) { // retune momentums
	    mom1[mom_ind] += x;
	    mom2[mom_ind] -= x;
            outgoingParticles[tune_particles.first ].setMomentum(mom1);
            outgoingParticles[tune_particles.second].setMomentum(mom2);
            out_mom = getTotalOutputMomentum();
	    std::sort(outgoingParticles.begin(), outgoingParticles.end(), G4ParticleLargerEkin());
            mon_non_cons = ini_mom - out_mom;
            pnc = mon_non_cons.rho();
	    if(verboseLevel > 2){
	      G4cout << " momentum non conservation tuning: " << G4endl 
		     << " e " << mon_non_cons.e() << " p " << pnc << G4endl;
	    }
            if(std::fabs(mon_non_cons.e()) < accuracy || pnc < accuracy) on_shell = true;
	  }
	  else {
	    if(verboseLevel > 2){
	      G4cout << " no appropriate solution found " << G4endl;
	    }
	  }; 
        }
	else {
	  if(verboseLevel > 2){
	    G4cout << " DET < 0 " << G4endl;
	  }
        };   
      }
      else { 
	if(verboseLevel > 2){
	  G4cout << " tuning impossible " << G4endl;
	}
      };   
    }
    else {
      on_shell = true;
    };
  }
  else {
    on_shell = true;
  };
  
}

std::pair<std::pair<G4int, G4int>, G4int> G4CollisionOutput::selectPairToTune(G4double de) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4CollisionOutput::selectPairToTune" << G4endl;
  }

  std::pair<G4int, G4int> tup(-1, -1);
  G4int i3 = -1; 
  std::pair<std::pair<G4int, G4int>, G4int> badp(tup, i3);

  if(outgoingParticles.size() < 2) {

    return badp;
  }
  else {

    G4int ibest1 = -1;
    G4int ibest2 = -1;  
    G4double pbest = 0.0;
    G4double pcut = 0.3 * std::sqrt(1.88 * std::fabs(de));
    G4double p1 = 0.0;
    G4double p2;
   
    for(G4int i = 0; i < G4int(outgoingParticles.size()) - 1; i++) {
      G4LorentzVector mom1 = outgoingParticles[i].getMomentum();

      for(G4int j = i+1; j < G4int(outgoingParticles.size()); j++) {
	G4LorentzVector mom2 = outgoingParticles[j].getMomentum();

	for(G4int l = G4LorentzVector::X; l<=G4LorentzVector::Z; l++) {
	  if(mom1[l] * mom2[l] < 0.0) { 
	    if(std::fabs(mom1[l]) > pcut && std::fabs(mom2[l]) > pcut) {

	      G4double psum = std::fabs(mom1[l]) + std::fabs(mom2[l]);  

	      if(psum > pbest) {
		ibest1 = i;
		ibest2 = j;
		i3 = l;
		p1 = mom1[l];
		p2 = mom2[l];
		pbest = psum;
	      };
	    };
	  };
	};
      };
    };    	       
    if(i3 > 0) {
      if(de > 0.0) {
	if(p1 > 0.0) {
	  tup.first  = ibest1;
	  tup.second = ibest2;
	}
        else {
	  tup.first  = ibest2;
	  tup.second = ibest1;
	};		 
      }
      else {
	if(p1 < 0.0) {
	  tup.first  = ibest2;
	  tup.second = ibest1;
	}
        else {
	  tup.first  = ibest1;
	  tup.second = ibest2;
	};		 
      }; 

      return std::pair<std::pair<G4int, G4int>, G4int>(tup, i3);
    };
  };  

  return badp;
}
