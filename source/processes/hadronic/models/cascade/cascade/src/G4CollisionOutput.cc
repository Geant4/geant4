#include "G4CollisionOutput.hh"
#include "G4ParticleLargerEkin.hh"
#include "algorithm"

typedef G4std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4CollisionOutput::G4CollisionOutput()
  : verboseLevel(2) {
  
  if (verboseLevel > 3) {
    G4cout << " >>> G4CollisionOutput::G4CollisionOutput" << G4endl;
  }
}

void G4CollisionOutput::setOnShell(G4InuclParticle* bullet, 
				   G4InuclParticle* target) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4CollisionOutput::setOnShell" << G4endl;
  }

  const G4double accuracy = 0.00001; // momentum concerves at the level of 1 eV

  on_shell = false;
    
  G4std::vector<G4double> ini_mom = bullet->getMomentum();
  G4std::vector<G4double> momt = target->getMomentum();

  for(G4int i = 0; i < 4; i++) ini_mom[i] += momt[i];
  
  G4std::vector<G4double> out_mom = getTotalOutputMomentum();
  G4std::vector<G4double> mon_non_cons(4);

  for(i = 0; i < 4; i++) mon_non_cons[i] = ini_mom[i] - out_mom[i];

  G4double pnc = sqrt(mon_non_cons[1] * mon_non_cons[1] + 
		      mon_non_cons[2] * mon_non_cons[2] + 
		      mon_non_cons[3] * mon_non_cons[3]);

  setRemainingExitationEnergy();       

  if(verboseLevel > 2){
    printCollisionOutput();
    G4cout << " momentum non conservation: " << endl << " e " << mon_non_cons[0]
	   << " p " << pnc << G4endl;
    G4cout << " remaining exitation " << eex_rest << G4endl;
  }

  if(fabs(mon_non_cons[0]) > accuracy || pnc > accuracy) { // renormalization
    G4int npart = outgoingParticles.size();
 
    if(npart > 0) {

      G4std::vector<G4double> last_mom = outgoingParticles[npart - 1].getMomentum(); 

      for(G4int i = 1; i < 4; i++) last_mom[i] += mon_non_cons[i];
      outgoingParticles[npart - 1].setMomentum(last_mom);
    }
    else {

      G4int nnuc = nucleiFragments.size();
 
      if(nnuc > 0) {

        G4std::vector<G4double> last_mom = nucleiFragments[nnuc - 1].getMomentum();

        for(G4int i = 1; i < 4; i++) last_mom[i] += mon_non_cons[i];
	nucleiFragments[nnuc - 1].setMomentum(last_mom);
	nucleiFragments[nnuc - 1].setEnergy();
      };
    }; 
    out_mom = getTotalOutputMomentum();
    for(G4int i = 0; i < 4; i++) mon_non_cons[i] = ini_mom[i] - out_mom[i];
    pnc = sqrt(mon_non_cons[1] * mon_non_cons[1] + 
	       mon_non_cons[2] * mon_non_cons[2] + 
	       mon_non_cons[3] * mon_non_cons[3]);

    if(verboseLevel > 2){
      printCollisionOutput();
      G4cout << " momentum non conservation after (1): " << G4endl 
	     << " e " << mon_non_cons[0] << " p " << pnc << G4endl;
    }
    G4bool need_hard_tuning = true;
    
    if(nucleiFragments.size() > 0) {
      for(G4int i = 0; i < G4int(nucleiFragments.size()); i++) {

        G4double eex = nucleiFragments[i].getExitationEnergyInGeV();

	if(eex > 0.0 && eex + mon_non_cons[0] >= 0.0) {
	  nucleiFragments[i].setExitationEnergy(1000.0 * (eex + mon_non_cons[0]));
	  need_hard_tuning = false;
	  break;
	};
      };
      if(need_hard_tuning && mon_non_cons[0] > 0.) {
	nucleiFragments[0].setExitationEnergy(1000.0 * mon_non_cons[0]);
	need_hard_tuning = false;
      };
    };
    
    if(need_hard_tuning) {
     
      pair<pair<G4int, G4int>, G4int> tune_par = selectPairToTune(mon_non_cons[0]);
      pair<G4int, G4int> tune_particles = tune_par.first;
      G4int mom_ind = tune_par.second;

      if(verboseLevel > 2) {
	G4cout << " p1 " << tune_particles.first << " p2 " << tune_particles.second
	       << " ind " << mom_ind << G4endl;
      }
      if(tune_particles.first >= 0 && tune_particles.second >= 0 &&
	 mom_ind >= 1) { // tunning possible

        G4std::vector<G4double> mom1 = outgoingParticles[tune_particles.first].getMomentum();
        G4std::vector<G4double> mom2 = outgoingParticles[tune_particles.second].getMomentum();
        G4double newE12 = mom1[0] + mom2[0] + mon_non_cons[0];
        G4double R = 0.5 * (newE12 * newE12 + mom2[0] * mom2[0] - mom1[0] * mom1[0]) / newE12;
        G4double Q = -(mom1[mom_ind] + mom2[mom_ind]) / newE12;
        G4double UDQ = 1.0 / (Q * Q - 1.0);
        G4double W = (R * Q + mom2[mom_ind]) * UDQ;
        G4double V = (mom2[0] * mom2[0] - R * R) * UDQ;
        G4double DET = W * W + V;

        if(DET > 0.0) {

          G4double x1 = -(W + sqrt(DET));
   	  G4double x2 = -(W - sqrt(DET));
	  // choose the appropriate solution
          G4bool xset = false;
 	  G4double x = 0.0;

	  if(mon_non_cons[0] > 0.0) { // x has to be > 0.0
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
	    sort(outgoingParticles.begin(), outgoingParticles.end(), G4ParticleLargerEkin());
            for(G4int i = 0; i < 4; i++) mon_non_cons[i] = ini_mom[i] - out_mom[i];
            pnc = sqrt(mon_non_cons[1] * mon_non_cons[1] + 
		       mon_non_cons[2] * mon_non_cons[2] + 
		       mon_non_cons[3] * mon_non_cons[3]);
	    if(verboseLevel > 2){
	      G4cout << " momentum non conservation tuning: " << G4endl 
		     << " e " << mon_non_cons[0] << " p " << pnc << G4endl;
	    }
            if(fabs(mon_non_cons[0]) < accuracy || pnc < accuracy) on_shell = true;
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

pair<pair<G4int, G4int>, G4int> G4CollisionOutput::selectPairToTune(G4double de) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4CollisionOutput::selectPairToTune" << G4endl;
  }

  pair<G4int, G4int> tup(-1, -1);
  G4int i3 = -1; 
  pair<pair<G4int, G4int>, G4int> badp(tup, i3);

  if(outgoingParticles.size() < 2) {

    return badp;
  }
  else {

    G4int ibest1 = -1;
    G4int ibest2 = -1;  
    G4double pbest = 0.0;
    G4double pcut = 0.3 * sqrt(1.88 * fabs(de));
    G4double p1 = 0.0;
    G4double p2;
   
    for(G4int i = 0; i < G4int(outgoingParticles.size()) - 1; i++) {

      G4std::vector<G4double> mom1 = outgoingParticles[i].getMomentum();

      for(G4int j = i+1; j < G4int(outgoingParticles.size()); j++) {

	G4std::vector<G4double> mom2 = outgoingParticles[j].getMomentum();

	for(G4int l = 1; l < 4; l++) {
	  if(mom1[l] * mom2[l] < 0.0) { 
	    if(fabs(mom1[l]) > pcut && fabs(mom2[l]) > pcut) {

	      G4double psum = fabs(mom1[l]) + fabs(mom2[l]);  

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

      return pair<pair<G4int, G4int>, G4int>(tup, i3);
    };
  };  

  return badp;
}
