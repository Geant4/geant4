#include "G4Collider.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4ParticleLargerEkin.hh"
#include "algorithm"

typedef vector<G4InuclElementaryParticle>::iterator particleIterator;

G4ElementaryParticleCollider::G4ElementaryParticleCollider()
  : verboseLevel(1) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::G4ElementaryParticleCollider" << G4endl;
  }
};

G4CollisionOutput  G4ElementaryParticleCollider::collide(G4InuclParticle* bullet,
							 G4InuclParticle* target) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::collide" << G4endl;
  }

  vector<G4double> totscm(4, 0.0); //::: fix
  vector<G4double> totlab(4, 0.0);


  //  generate nucleon or pion collission with NUCLEON 
  //  or pion with quasideutron

  if(verboseLevel > 2){
    G4cout << " here " << G4endl;
  }

  G4CollisionOutput output;  
  G4InuclElementaryParticle* particle1 =
    dynamic_cast<G4InuclElementaryParticle*>(bullet);
  G4InuclElementaryParticle* particle2 =	
    dynamic_cast<G4InuclElementaryParticle*>(target);
	     
  if(particle1 && particle2) { // particle / particle 
    if(!particle1->photon() && !particle2->photon()) { // ok
      if(particle1->nucleon() || particle2->nucleon()) { // ok

	if(verboseLevel > 2){
	  G4cout << " here1 " << G4endl;
	  particle1->printParticle();
	  particle2->printParticle();

	  vector<G4double> momb = particle1->getMomentum();
	  vector<G4double> momt = particle2->getMomentum();

	  for(G4int i = 0; i < 4; i++) momb[i] += momt[i];
	  G4cout << " total input: px " << momb[1] << " py " << momb[2] 
		 << " pz " << momb[3] << " e " << momb[0] << endl;
	}
	G4LorentzConvertor convertToSCM;

	if(particle2->nucleon()) {
          convertToSCM.setBullet(particle1->getMomentum(), particle1->getMass());
          convertToSCM.setTarget(particle2->getMomentum(), particle2->getMass());
	}
	else {
          convertToSCM.setBullet(particle2->getMomentum(), particle2->getMass());
          convertToSCM.setTarget(particle1->getMomentum(), particle1->getMass());
	};  
        convertToSCM.toTheCenterOfMass();
 
        G4double ekin = convertToSCM.getKinEnergyInTheTRS();
        G4double etot_scm = convertToSCM.getTotalSCMEnergy();
        G4double pscm = convertToSCM.getSCMMomentum();

	if(verboseLevel > 2){
	  G4cout << " ekin " << ekin << " etot_scm " << etot_scm << " pscm " <<
	    pscm << G4endl;
	}
        vector<G4InuclElementaryParticle> particles = 	    
	  generateSCMfinalState(ekin, etot_scm, pscm, particle1, particle2, &convertToSCM);

	if(verboseLevel > 2){
	  G4cout << " particles " << particles.size() << G4endl;
	  for(G4int i = 0; i < particles.size(); i++) 
	    particles[i].printParticle();
	}
	cout << "1" << endl;
	if(!particles.empty()) { // convert back to Lab
	  cout << "2" << endl;

	  /*	
		if(verboseLevel > 2){
		vector<G4double> totscm(4, 0.0); // moded as private variables
		cout << "3" << endl;
		vector<G4double> totlab(4, 0.0);
		cout << "4" << endl;
		}
	
	  */


	  particleIterator ipart;
	  cout << "5" << endl;
	  for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
	    cout << "6" << endl;
	
	    if(verboseLevel > 2){
	      cout << "7" << endl;
	      vector<G4double> mom_scm = ipart->getMomentum();
	      cout << "8" << endl;
	      cout << mom_scm[0] << " " <<  mom_scm[1] << " " <<  mom_scm[2] << " " <<  mom_scm[3] << endl;

	      for(G4int i = 0; i < 4; i++) {
		totscm[i] += mom_scm[i];
		cout << "8" << "/" << i << endl;
	      }
	    }
	    cout << "9" << endl;
	    vector<G4double> mom = 
	      convertToSCM.backToTheLab(ipart->getMomentum());
	    cout << "10" << endl;

	
	    if(verboseLevel > 2){
	      for(G4int i = 0; i < 4; i++) totlab[i] += mom[i];
	    }
       

	    ipart->setMomentum(mom); 
	  };

	  sort(particles.begin(), particles.end(), G4ParticleLargerEkin());

	  if(verboseLevel > 2){
	    G4cout << " In SCM: total outgoing momentum " << G4endl 
		   << " E " << totscm[0] << " px " << totscm[1]
		   << " py " << totscm[2] << " pz " << totscm[3] << G4endl; 
	    G4cout << " In Lab: total outgoing momentum " << G4endl 
		   << " E " << totlab[0] << " px " << totlab[1]
		   << " py " << totlab[2] << " pz " << totlab[3] << G4endl; 
	  }
	};
	
	output.addOutgoingParticles(particles);
      }
      else {
        if(particle1->quasi_deutron() || particle2->quasi_deutron()) {
	  if(particle1->pion() || particle2->pion()) {

	    G4LorentzConvertor convertToSCM;

            if(particle1->pion()) {
              convertToSCM.setBullet(particle1->getMomentum(), particle1->getMass());
              convertToSCM.setTarget(particle2->getMomentum(), particle2->getMass());
	    }
	    else {
              convertToSCM.setBullet(particle2->getMomentum(), particle2->getMass());
              convertToSCM.setTarget(particle1->getMomentum(), particle1->getMass());
	    }; 
            convertToSCM.toTheCenterOfMass(); 

            G4double etot_scm = convertToSCM.getTotalSCMEnergy();

	    if(verboseLevel > 2){
	      G4cout << " etot_scm " << etot_scm << G4endl;
	    }
            vector<G4InuclElementaryParticle> particles = 
	      generateSCMpionAbsorption(etot_scm, particle1, particle2);

	    if(verboseLevel > 2){
	      G4cout << " particles " << particles.size() << G4endl;

	      for(G4int i = 0; i < particles.size(); i++) 
		particles[i].printParticle();
	    }
	    if(!particles.empty()) { // convert back to Lab

	      /*
		if(verboseLevel > 2){
		vector<G4double> totscm(4, 0.0);
		}
	      */

	      particleIterator ipart;

	      for(ipart = particles.begin(); ipart != particles.end(); ipart++) {

		if(verboseLevel > 2){

		  vector<G4double> mom_scm = ipart->getMomentum();

		  for(G4int i = 0; i < 4; i++) totscm[i] += mom_scm[i];
		}

	        vector<G4double> mom = 
	          convertToSCM.backToTheLab(ipart->getMomentum());

	        ipart->setMomentum(mom); 
	      };
	      sort(particles.begin(), particles.end(), G4ParticleLargerEkin());

	      if(verboseLevel > 2){
		G4cout << " In SCM: total outgoing momentum " << G4endl 
		       << " E " << totscm[0] << " px " << totscm[1]
		       << " py " << totscm[2] << " pz " << totscm[2] << G4endl; 
	      }
	
 	      output.addOutgoingParticles(particles);
            };
	  }
	  else {
	    G4cout << " ElementaryParticleCollider -> can collide just pions with deutron at the moment " << G4endl;
	  }; 
	}
	else {
	  G4cout << " ElementaryParticleCollider -> can collide just smth. with nucleon or deutron at the moment " << G4endl;
        };
      };  
    }
    else {
      G4cout << " ElementaryParticleCollider -> can not collide photons at the moment " << G4endl;
    }; 
  }
  else {
    G4cout << " ElementaryParticleCollider -> can collide only particle with particle " << G4endl;
  }; 	 	 

  return output;
}

G4int G4ElementaryParticleCollider::generateMultiplicity(G4int is, 
							 G4double ekin) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::generateMultiplicity" << G4endl;
  }

  const G4double asig[4][6][31] = {
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 24.3, 24.1, 24.0, 26.3,
    28.6, 24.8, 19.9, 19.2, 17.4, 15.3, 13.5, 12.3, 11.9, 10.4,
    11.8, 11.4, 11.0, 10.8, 10.9, 11.7, 11.4, 10.2, 11.0, 11.0,
    9.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.45, 2.90,
    4.10, 5.30, 22.0, 21.2, 19.9, 14.9, 12.6, 12.3, 11.3, 10.8,
    9.50, 8.27, 7.20, 6.70, 6.25, 6.04, 5.89, 5.70, 5.60, 5.05,
    4.17, 4.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 2.27, 4.28, 7.81, 8.03, 8.11, 7.90, 7.82,
    7.61, 7.47, 7.07, 7.66, 7.05, 6.71, 6.38, 6.36, 6.37, 6.57,
    6.01, 5.48, 6.80, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 1.20, 2.85, 3.70, 4.81, 5.33,
    7.74, 6.91, 6.94, 7.57, 7.21, 7.11, 7.10, 6.93, 6.79, 6.71,
    6.55, 6.55, 6.15, 8.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.005,0.54, 0.74,
    0.86, 0.91, 1.10, 1.16, 1.36, 1.40, 1.43, 1.47, 1.47, 1.43,
    1.38, 1.38, 1.63, 1.36, 2.80, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 34.0, 46.2, 46.9, 45.2, 47.1,
    42.3, 41.8, 41.2, 41.6, 41.6, 41.0, 43.0, 42.4, 40.0, 39.9,
    39.8, 42.0, 40.0, 39.8, 39.6, 38.7,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 33.0, 31.3, 29.5, 27.8,
    14.6, 16.0, 17.5, 18.3, 19.4, 18.7, 15.6, 14.8, 13.6, 12.5,
    12.2, 11.9, 11.4, 11.2, 10.1, 9.62, 8.41, 7.14, 7.09, 5.04,
    10.2, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.00, 4.00,
    6.50, 19.5, 20.8, 19.8, 18.6, 17.7, 14.4, 13.5, 10.4, 10.1,
    12.0, 8.87, 8.51, 8.49, 9.20, 8.29, 7.43, 8.20, 4.69, 4.57,
    4.06, 4.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.68, 2.76, 3.85, 8.35, 12.7, 12.3, 10.8, 12.0,
    10.9, 10.2, 10.6, 12.2, 11.8, 11.0, 10.4, 10.5, 11.0, 10.6,
    11.8, 12.5, 13.2, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.22, 0.52, 0.64, 1.19, 1.52, 1.75, 1.51,
    2.04, 1.85, 1.70, 1.92, 1.66, 1.74, 1.50, 1.39, 1.35, 1.41,
    1.48, 1.43, 1.35, 3.20, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.14, 0.24, 0.30, 0.46,
    0.85, 1.40, 1.54, 1.52, 1.47, 1.48, 1.49, 1.42, 1.39, 1.37,
    1.22, 1.19, 0.93, 0.00, 2.10, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 35.0, 40.0, 42.4, 42.3, 41.0,
    40.9, 40.4, 39.8, 35.0, 33.6, 41.2, 41.0, 41.1, 41.2, 41.2,
    39.6, 36.0, 36.0, 36.2, 0.00, 40.2,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 60.0, 38.0, 30.6, 24.0,
    18.5, 12.8, 13.6, 9.15, 8.20, 7.80, 7.10, 6.40, 5.81, 5.85,
    5.50, 5.33, 5.40, 5.50, 4.90, 5.02, 5.00, 4.98, 4.96, 4.96,
    4.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.38, 0.75,
    1.28, 2.26, 11.8, 11.1, 7.56, 6.04, 5.68, 4.15, 3.21, 2.12,
    1.66, 1.54, 1.53, 1.47, 1.33, 1.39, 1.35, 1.29, 1.23, 1.13,
    1.06, 0.70, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.05, 2.19, 6.02, 7.36, 6.97, 5.83, 6.53, 5.09,
    4.24, 3.24, 3.31, 3.11, 2.91, 2.78, 2.72, 2.68, 2.48, 2.27,
    2.02, 1.77, 4.40, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.06, 0.35, 0.55, 0.77, 0.89, 0.80, 0.86,
    0.97, 0.90, 0.86, 0.84, 0.84, 0.83, 0.81, 0.79, 0.82, 0.76,
    0.88, 0.89, 0.89, 3.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.02, 0.07, 0.33, 0.92, 1.39,
    2.11, 1.81, 2.39, 2.60, 2.19, 1.70, 1.60, 0.68, 1.43, 1.46,
    1.46, 1.37, 1.16, 1.09, 2.60, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 18.9, 27.2, 34.9, 29.1, 30.8,
    29.6, 28.2, 27.5, 26.9, 26.3, 25.9, 25.6, 25.2, 26.1, 25.5,
    25.4, 25.3, 25.1, 24.9, 24.8, 24.1,
    5.90, 9.40, 24.5, 62.6, 65.3, 41.3, 29.3, 24.3, 22.7, 22.9,
    23.2, 28.4, 11.7, 10.1, 8.30, 7.16, 6.49, 6.36, 6.60, 5.84,
    5.30, 4.50, 3.90, 4.40, 4.74, 0.794,0.824,0.714,0.59, 0.00,
    4.60, 0.00, 0.00, 0.00, 0.00, 0.10, 0.40, 2.70, 3.50, 5.30,
    6.60, 9.10, 17.6, 12.2, 9.78, 7.51, 6.91, 6.86, 6.46, 6.19,
    5.13, 3.90, 2.82, 3.10, 3.12, 2.52, 2.22, 2.02, 2.01, 1.98,
    2.14, 1.20, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.76, 2.63, 3.72, 6.53, 7.47, 7.94, 7.12, 6.85,
    6.09, 5.35, 4.12, 3.85, 3.68, 4.09, 3.58, 3.29, 3.08, 2.93,
    2.80, 2.65, 3.30, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.59, 0.74, 1.47, 4.10, 4.78, 4.90, 5.07,
    5.50, 5.48, 5.03, 4.65, 4.39, 4.06, 3.53, 3.08, 3.05, 2.91,
    3.42, 3.93, 3.93, 4.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.01,0.007, 0.03,0.099,0.251,.0376,
    0.419, 0.582, 0.755, 0.777, 1.13, 1.08, 1.13, 1.08, 0.962, 0.866,
    0.738, 0.674, 0.645, 0.613, 1.30, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 31.3, 46.0, 30.0, 35.7, 33.4,
    31.6, 30.4, 29.6, 28.9, 28.5, 28.1, 27.5, 31.0, 27.7, 27.8,
    26.1, 25.2, 6.92, 6.70, 0.00, 25.7
  };

  const G4double large_cut = 4.0;
  pair<G4int, G4double> iksk = getPositionInEnergyScale2(ekin);
  G4int ik = iksk.first;
  G4double sk = iksk.second;
  G4int l = is;

  if(l == 4) l = 1; 
  if(l == 10) l = 3; 
  if(l == 5 || l == 6) l = 4; 

  vector<G4double> sigm(5);
  G4double stot = 0.0;

  if(l == 7 || l == 14) { // pi0 P or pi0 N
    for(G4int j = 0; j < 5; j++) {
      sigm[j] = fabs(0.5 * (asig[2][j][ik - 1] + asig[3][j][ik - 1] +
			    sk * (asig[2][j][ik] + asig[3][j][ik] - 
				  asig[2][j][ik - 1] - asig[3][j][ik - 1])));
      stot += sigm[j];
    };
  }
  else {
    for(G4int j = 0; j < 5; j++) {
      sigm[j] = fabs(asig[l - 1][j][ik - 1] + sk * (asig[l - 1][j][ik] 
						    - asig[l - 1][j][ik - 1]));
      stot += sigm[j];
    };
  };
  G4double sl = inuclRndm();
  G4double ptot = 0.0;
  G4int mul;

  for(G4int i = 0; i < 5; i++) {
    ptot += sigm[i] / stot;
    if(sl <= ptot) {
      mul = i;

      break;
    };  
  };
  if(ekin > large_cut && mul == 1) mul = 2;

  if(verboseLevel > 3){
    G4cout << " multiplicity " << mul + 2 << G4endl; 
  }

  return mul + 2;
}

vector<G4InuclElementaryParticle> G4ElementaryParticleCollider:: 
generateSCMfinalState(G4double ekin, 
		      G4double etot_scm, 
		      G4double pscm,
		      G4InuclElementaryParticle* particle1,
		      G4InuclElementaryParticle* particle2, 
		      G4LorentzConvertor* toSCM) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::generateSCMfinalState" << G4endl;
  }

  const G4double ang_cut = 0.9999;
  const G4double difr_const = 0.3678794;   
  const G4int itry_max = 10;
  G4InuclElementaryParticle dummy;
  vector<G4InuclElementaryParticle> particles;
  vector<G4int> particle_kinds;
  G4int type1 = particle1->type();
  G4int type2 = particle2->type();
  G4int is = type1 * type2;

  if(verboseLevel > 3){
    G4cout << " is " << is << G4endl;
  }

  G4int multiplicity = 0;
  G4bool generate = true;
   
  while(generate) {      
    if(multiplicity == 0) {
      multiplicity = generateMultiplicity(is, ekin);
    }
    else {
      multiplicity = generateMultiplicity(is, ekin);
      particle_kinds.resize(0);
    };
    if(multiplicity == 2) { // 2 -> 2

      G4int kw;

      if(reChargering(ekin, is)) { // rechargering
	kw = 2;
	switch (is) {
	case 6: // pi+ N -> pi0 P
	  particle_kinds.push_back(7);
	  particle_kinds.push_back(1);
	  break;    
	case 5: // pi- P -> pi0 N
	  particle_kinds.push_back(7);
	  particle_kinds.push_back(2);
	  break;    
	case 7: // pi0 P -> pi+ N
	  particle_kinds.push_back(3);
	  particle_kinds.push_back(2);
	  break;    
	case 14: // pi0 N -> pi- P
	  particle_kinds.push_back(5);
	  particle_kinds.push_back(1);
	  break;    
	default: 

	  G4cout << " strange recharge: " << is << G4endl;

	  particle_kinds.push_back(type1);
	  particle_kinds.push_back(type2);
	  kw = 1;
	};
      }
      else { // just elastic
	kw = 1;
	particle_kinds.push_back(type1);
	particle_kinds.push_back(type2);       
      };

      vector<G4double> mom;

      if(kw == 2) { // need to rescale momentum

	G4double m1 = dummy.getParticleMass(particle_kinds[0]);

	m1 *= m1;

	G4double m2 = dummy.getParticleMass(particle_kinds[1]);

	m2 *= m2;	 

	G4double a = 0.5 * (etot_scm * etot_scm - m1 - m2);
	G4double np = sqrt((a * a - m1 * m2) / (m1 + m2 + 2.0 * a));

	mom = particleSCMmomentumFor2to2(is, kw, ekin, np);
      }
      else {
	mom = particleSCMmomentumFor2to2(is, kw, ekin, pscm);
      };

      if(verboseLevel > 3){
	G4cout << " before rotation px " << mom[1] << " py " << mom[2] <<
	  " pz " << mom[3] << G4endl;
      }

      mom = toSCM->rotate(mom); 

      if(verboseLevel > 3){
	G4cout << " after rotation px " << mom[1] << " py " << mom[2] <<
	  " pz " << mom[3] << G4endl;
      }
      vector<G4double> mom1(4);

      for(G4int i = 1; i < 4; i++) mom1[i] = -mom[i];
      particles.push_back(G4InuclElementaryParticle(mom, particle_kinds[0]));
      particles.push_back(G4InuclElementaryParticle(mom1, particle_kinds[1]));
      generate = false;
    } else { // 2 -> many
      particle_kinds = generateOutgoingKindsFor2toMany(is, multiplicity, ekin);

      G4int itry = 0;
      G4bool bad = true;
      G4int knd_last = particle_kinds[multiplicity - 1];
      G4double mass_last = dummy.getParticleMass(knd_last);

      if(verboseLevel > 3){
	G4cout << " knd_last " << knd_last << " mass " << mass_last << G4endl;
      }

      while(bad && itry < itry_max) {
	itry++;

	if(verboseLevel > 3){
	  G4cout << " itry in while " << itry << G4endl;
	}

	vector<G4double> modules = 
	  generateMomModules(particle_kinds, multiplicity, is, ekin, etot_scm);
	if(modules.size() == multiplicity) {
	  if(multiplicity == 3) { 
	    vector<G4double> mom3 = 
	      particleSCMmomentumFor2to3(is, knd_last, ekin, modules[2]);
	    mom3 = toSCM->rotate(mom3);
	    // generate the momentum of first
	    G4double ct = -0.5 * (modules[2] * modules[2] + 
				  modules[0] * modules[0] - 
				  modules[1] * modules[1]) /
	      modules[2] / modules[0];   
	    if(fabs(ct) < ang_cut) {

	      if(verboseLevel > 2){
		G4cout << " ok for mult " << multiplicity << G4endl;
	      }

	      vector<G4double> mom1 = generateWithFixedTheta(ct, modules[0]);

	      mom1 = toSCM->rotate(mom3, mom1);

	      vector<G4double> mom2(4);

	      for(G4int i = 1; i < 4; i++) mom2[i] = - (mom3[i] + mom1[i]);
	      bad = false;
	      generate = false;
	      particles.push_back(G4InuclElementaryParticle(mom1, particle_kinds[0]));
	      particles.push_back(G4InuclElementaryParticle(mom2, particle_kinds[1]));
	      particles.push_back(G4InuclElementaryParticle(mom3, particle_kinds[2]));
	    };
	  }
	  else { // multiplicity > 3
	    //        generate first mult - 2 momentums
	    vector<vector<G4double> > scm_momentums;
	    vector<G4double> tot_mom(4);

	    for(G4int i = 0; i < multiplicity - 2; i++) {

	      G4double p0 = particle_kinds[i] < 3 ? 0.36 : 0.25;
	      G4double alf = 1.0 / p0 / (p0 - (modules[i] + p0) *
					 exp(-modules[i] / p0));
	      G4double st = 2.0;
	      G4int itry1 = 0;

	      while(fabs(st) > ang_cut && itry1 < itry_max) {
		itry1++;

		G4double s1 = modules[i] * inuclRndm();
		G4double s2 = alf * difr_const * p0 * inuclRndm();

		if(verboseLevel > 3){
		  G4cout << " s1 * alf * exp(-s1 / p0) " << s1 * alf * exp(-s1 / p0) 
			 << " s2 " << s2 << G4endl;
		}

		if(s1 * alf * exp(-s1 / p0) > s2) st = s1 / modules[i];
	      }; 

	      if(verboseLevel > 3){
		G4cout << " itry1 " << itry1 << " i " << i << " st " << st << G4endl;
	      }
	      if(itry1 == itry_max) {

		if(verboseLevel > 2){
		  G4cout << " high energy angles generation: itry1 " << itry1 << G4endl;
		}

		st = 0.5 * inuclRndm();
	      };

	      G4double ct = sqrt(1.0 - st * st);

	      if(inuclRndm() > 0.5) ct = -ct;

	      G4double pt = modules[i]*st;
	      G4double phi = randomPHI();
	      vector<G4double> mom(4);

	      mom[1] = pt * cos(phi);
	      mom[2] = pt * sin(phi);
	      mom[3] = modules[i] * ct;
	      for(G4int i = 1; i < 4; i++) tot_mom[i] += mom[i];		 
	      scm_momentums.push_back(mom);
	    }; 
	    //         handle last two
	    G4double tot_mod = sqrt(tot_mom[1] * tot_mom[1] + 
				    tot_mom[2] * tot_mom[2] + tot_mom[3] * tot_mom[3]); 
	    G4double ct = -0.5 * (tot_mod * tot_mod + 
				  modules[multiplicity - 2] * modules[multiplicity - 2] -
				  modules[multiplicity - 1] * modules[multiplicity - 1]) / tot_mod /
	      modules[multiplicity - 2];  

	    if(verboseLevel > 2){
	      G4cout << " ct last " << ct << G4endl;
	    }            

	    if(fabs(ct) < ang_cut) {
	      for(G4int i = 0; i < multiplicity - 2; i++) 
		scm_momentums[i] = toSCM->rotate(scm_momentums[i]);
	      tot_mom = toSCM->rotate(tot_mom);  

	      vector<G4double> mom = 
		generateWithFixedTheta(ct, modules[multiplicity - 2]);

	      mom = toSCM->rotate(tot_mom, mom);
	      scm_momentums.push_back(mom);
	      // and the last one
	      vector<G4double> mom1(4);

	      for(i = 1; i < 4; i++) mom1[i] = -mom[i] - tot_mom[i];
	      scm_momentums.push_back(mom1);  
	      bad = false;
	      generate = false;

	      if(verboseLevel > 2){
		G4cout << " ok for mult " << multiplicity << G4endl;
	      }

	      for(i = 0; i < multiplicity; i++) {
		particles.push_back(G4InuclElementaryParticle(
							      scm_momentums[i], particle_kinds[i]));
	      };
	    };
	  }; 
	};
      };
      if(itry == itry_max) {

	if(verboseLevel > 2){
	  G4cout << " cannot generate the distr. for mult " << multiplicity  <<
	    G4endl << " and set it to " << multiplicity - 1 << G4endl;
	}

      };
    };
  }; 

  if (verboseLevel > 3) {
    G4cout << " <<< G4ElementaryParticleCollider::generateSCMfinalState" << G4endl;
  }

  return particles;
}

vector<G4double> G4ElementaryParticleCollider::
generateMomModules(
		   const vector<G4int>& kinds, 
		   G4int mult, 
		   G4int is, 
		   G4double ekin, 
		   G4double etot_cm) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::generateMomModules" << G4endl;
  }

  if(verboseLevel > 2){
    G4cout << " mult " << mult << " is " << is << " ekin " << ekin << " etot_cm " <<
      etot_cm << G4endl;
  }

  const G4int itry_max = 10;
  const G4double small = 1.e-10;
  G4InuclElementaryParticle dummy;
  G4int itry = 0;
  G4double sr;
  vector<G4double> modules(mult);
  vector<G4double> masses2(mult);

  for(G4int i = 0; i < mult; i++) {
    G4double mass = dummy.getParticleMass(kinds[i]);
    masses2[i] = mass * mass;
  };

  G4double mass_last = sqrt(masses2[mult - 1]);

  if(verboseLevel > 3){
    G4cout << " knd_last " << kinds[mult - 1] << " mlast " << mass_last << G4endl;
  }

  while (itry < itry_max) {
    itry++;
    if(verboseLevel > 3){
      G4cout << " itry in generateMomModules " << itry << G4endl;
    }

    G4int ilast = -1;
    G4double eleft = etot_cm;

    for(G4int i = 0; i < mult - 1; i++) {

      G4double pmod = 
	getMomModuleFor2toMany(is, mult, kinds[i], ekin);

      if(pmod < small) break;
      eleft -= sqrt(pmod * pmod + masses2[i]);

      if(verboseLevel > 3){
	G4cout << " kp " << kinds[i] << " pmod " << pmod << " mass2 " << masses2[i] << G4endl;
	G4cout << " x1 " << eleft - mass_last << G4endl;
      }

      if(eleft <= mass_last) break;
      ilast++;
      modules[i] = pmod;
    };

    if(ilast == mult - 2) {

      G4double plast = eleft * eleft - masses2[mult - 1];

      if(verboseLevel > 2){
	G4cout << " plast ** 2 " << plast << G4endl;
      }

      if(plast > small) {
	plast = sqrt(plast);
	modules[mult - 1] = plast;      
	if(mult == 3) { 
	  if(satisfyTriangle(modules)) {

	    return modules;
	  }
	}
	else {

	  return modules;
	}; 	 
      };
    };
  };

  modules.resize(0);

  return modules;    
}

G4bool G4ElementaryParticleCollider::satisfyTriangle(
						     const vector<G4double>& modules) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::satisfyTriangle" << G4endl;
  }

  G4bool good = true;

  if(modules.size() == 3) {
    if(fabs(modules[1] - modules[2]) > modules[0] || 
       modules[0] > modules[1] + modules[2] ||
       fabs(modules[0] - modules[2]) > modules[1] ||
       modules[1] > modules[0] + modules[2] ||
       fabs(modules[0] - modules[1]) > modules[2] ||
       modules[2] > modules[1] + modules[0]) good = false;
  };

  return good;
}
      
G4int G4ElementaryParticleCollider::getIL(G4int is, 
					  G4int mult) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::getIL" << G4endl;
  }

  const G4int ifdef[4][7] = {
    {2, 3, 2, 2, 3, 3, 2}, 
    {4, 4, 3, 4, 4, 4, 3}, 
    {5, 6, 4, 5, 5, 5, 4},
    {7, 7, 5, 7, 6, 6, 5}
  };

  G4int l = is;

  if(l == 14) {
    l = 5;
  }
  else if(l == 7) {
    l = 6;
  }
  else if(l == 10) {
    l = 7;
  };

  return ifdef[mult - 3][l - 1];
}

vector<G4int> G4ElementaryParticleCollider::
generateOutgoingKindsFor2toMany(
				G4int is, 
				G4int mult, 
				G4double ekin) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::generateOutgoingKindsFor2toMany" << G4endl;
  }

  const G4double bsig[4][20][20] = {
    1.20,3.70,3.98,3.85,3.51,2.90,2.86,2.81,2.77,2.80,
    2.54,2.00,1.90,1.75,1.68,1.61,1.54,1.40,1.25,1.17,
    8.00,14.0,13.0,12.0,11.4,9.70,9.41,8.52,8.03,6.70,
    5.73,5.20,4.80,4.54,4.36,4.28,4.16,4.10,3.80,3.00,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.37,0.41,0.92,1.10,0.98,0.85,0.74,0.74,0.83,
    0.71,0.68,0.53,0.48,0.41,0.41,0.41,0.37,0.37,0.37,
    0.  ,0.30,1.22,2.51,2.67,2.95,2.95,2.95,2.96,2.84,
    2.80,2.70,3.00,2.80,2.46,2.46,2.45,2.40,2.40,2.20,
    0.  ,1.40,2.37,4.07,3.90,3.80,3.70,3.70,3.60,3.40,
    3.20,3.50,3.17,3.00,2.80,2.70,2.60,2.60,2.40,2.20,
    0.  ,0.20,0.28,0.31,0.36,0.38,0.40,0.42,0.43,0.44,
    0.46,0.48,0.50,0.80,0.71,0.80,0.96,1.20,1.00,0.91,
    0.  ,0.  ,0.14,0.14,0.14,0.14,0.12,0.11,0.11,0.11,
    0.14,0.17,0.16,0.16,0.16,0.15,0.15,0.15,0.15,0.14,
    0.  ,0.  ,0.02,0.21,0.74,1.10,1.50,1.76,1.98,2.40,
    2.50,2.60,2.60,2.50,2.50,2.40,2.40,2.30,2.10,2.00,
    0.  ,0.  ,0.10,0.40,1.15,1.60,1.80,2.19,2.80,2.30,
    2.90,2.60,2.60,2.60,2.50,2.50,2.50,2.40,2.20,1.85,
    0.  ,0.  ,0.80,1.90,1.50,1.80,1.70,1.60,1.80,1.90,
    1.80,1.60,1.50,1.43,1.41,1.30,1.28,1.30,1.60,1.70,
    0.  ,0.  ,0.14,0.16,0.17,0.17,0.21,0.21,0.22,0.23,
    0.24,0.25,0.41,0.36,0.40,0.46,0.52,0.50,0.50,0.46,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.04,0.06,0.05,0.04,0.04,0.04,
    0.04,0.05,0.05,0.05,0.04,0.04,0.04,0.03,0.03,0.03,
    0.  ,0.  ,0.  ,.005,0.02,0.07,0.11,0.15,0.17,0.19,
    0.24,0.25,0.26,0.26,0.25,0.25,0.24,0.24,0.23,0.21,
    0.  ,0.  ,0.  ,0.  ,0.05,0.09,0.18,0.20,0.22,0.28,
    0.36,0.41,0.42,0.43,0.44,0.45,0.46,0.46,0.47,0.47,
    0.  ,0.  ,0.  ,0.  ,0.09,0.16,0.18,0.22,0.28,0.23,
    0.29,0.26,0.26,0.26,0.25,0.25,0.25,0.24,0.22,0.18,
    0.  ,0.  ,0.  ,0.  ,0.15,0.18,0.17,0.16,0.18,0.19,
    0.18,0.16,0.15,0.14,0.14,0.13,0.12,0.13,0.16,0.17,
    0.  ,0.  ,0.  ,0.  ,0.02,0.02,0.02,0.02,0.04,0.05,
    0.05,0.06,0.07,0.10,0.11,0.10,0.09,0.09,0.09,0.09,
    0.  ,0.  ,0.  ,0.  ,0.17,0.16,0.15,0.16,0.17,0.18,
    0.20,0.21,0.22,0.23,0.24,0.21,0.18,0.19,0.20,0.21,
    12.0,14.0,13.1,12.0,11.0,10.0,8.70,6.30,5.20,7.80,
    6.40,6.20,6.20,6.80,5.40,4.80,4.80,2.60,2.50,1.90,
    1.90,2.50,2.60,1.80,1.70,1.68,1.61,1.50,1.46,1.40,
    1.31,1.17,1.14,1.10,1.09,1.03,1.00,0.99,0.97,0.96,
    5.60,4.30,4.20,4.80,5.00,2.70,3.20,2.60,3.40,2.80,
    1.16,1.14,1.15,1.30,1.80,1.60,1.20,1.10,1.10,1.20,
    0.12,1.40,1.30,1.20,2.40,2.00,1.60,3.80,3.40,3.20,
    3.50,3.80,3.70,3.30,3.30,3.40,3.50,3.16,3.80,4.10,
    0.  ,0.77,1.75,5.28,6.30,5.90,5.30,4.90,4.80,4.70,
    4.20,3.72,3.60,3.40,3.20,3.50,3.90,4.10,4.30,4.50,
    0.  ,0.16,0.35,0.91,2.80,3.60,3.20,2.60,1.90,1.14,
    1.15,1.16,1.14,1.11,1.08,0.99,0.94,0.91,0.86,0.83,
    0.56,0.43,0.42,0.96,1.20,0.80,0.64,0.71,0.84,1.20,
    1.80,3.50,3.40,3.20,2.80,2.60,2.70,2.40,2.80,3.10,
    0.01,0.02,0.01,0.01,0.02,0.10,0.16,0.38,0.34,0.32,
    0.35,0.38,0.37,0.33,0.35,0.36,0.31,0.42,0.41,0.39,
    0.  ,0.07,0.17,0.53,0.63,0.59,0.53,0.49,0.48,0.47,
    0.42,0.37,0.36,0.34,0.32,0.35,0.39,0.41,0.43,0.45,
    0.09,0.18,0.21,0.28,0.36,0.37,0.20,0.64,0.68,0.56,
    0.61,0.48,0.34,0.36,0.36,0.32,0.41,0.39,0.39,0.36,
    0.,0.02,0.04,0.18,0.28,0.36,0.32,0.26,0.19,0.12,
    0.12,0.13,0.14,0.15,0.12,0.11,0.10,0.09,0.07,0.04,
    0.10,0.12,0.11,0.10,0.12,0.16,0.14,0.13,0.14,0.12,
    0.15,0.18,0.21,0.21,0.14,0.13,0.12,0.11,0.09,0.08,
    0.02,0.11,0.10,0.09,0.11,0.17,0.16,0.14,0.12,0.11,
    0.13,0.12,0.11,0.11,0.10,0.09,0.08,0.06,0.04,0.03,
    0.  ,0.  ,0.04,0.07,0.08,0.11,0.12,0.16,0.15,0.15,
    0.14,0.13,0.12,0.18,0.17,0.16,0.15,0.09,0.09,0.03,
    0.  ,0.  ,0.09,0.10,0.11,0.12,0.28,0.49,0.58,0.53,
    0.48,0.46,0.46,0.44,0.43,0.43,0.42,0.39,0.39,0.31,
    0.  ,0.  ,0.01,0.05,0.06,0.09,0.11,0.12,0.12,0.14,
    0.15,0.16,0.15,0.14,0.14,0.13,0.18,0.16,0.14,0.13,
    0.  ,0.  ,0.  ,0.02,0.04,0.09,0.20,0.39,0.41,0.42,
    0.42,0.44,0.45,0.43,0.41,0.43,0.42,0.39,0.39,0.30,
    0.  ,0.  ,0.  ,0.  ,0.01,0.02,0.08,0.04,0.08,0.09,
    0.10,0.11,0.12,0.11,0.09,0.09,0.06,0.06,0.06,0.06,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.01,0.02,0.08,0.09,0.09,
    0.09,0.09,0.10,0.11,0.11,0.08,0.07,0.07,0.07,0.06,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.02,0.04,0.12,0.11,0.10,
    0.09,0.08,0.08,0.08,0.07,0.07,0.07,0.06,0.05,0.04,
    1.80,9.19,7.40,5.29,3.64,3.48,2.47,2.21,1.30,0.95,
    0.82,0.79,0.76,0.75,0.74,0.74,0.72,0.70,0.65,0.61,
    0.46,2.60,3.70,2.27,2.40,2.00,1.68,1.00,0.82,0.71,
    0.72,0.74,0.71,0.68,0.65,0.61,0.57,0.53,0.48,0.45,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.02,0.15,1.55,2.31,1.98,1.52,1.30,1.10,0.68,0.51,
    0.42,0.39,0.35,0.34,0.33,0.32,0.31,0.30,0.30,0.29,
    0.01,0.95,1.09,1.63,1.50,1.40,1.38,1.20,1.00,0.81,
    0.63,0.52,0.38,0.34,0.32,0.31,0.24,0.23,0.22,0.21,
    0.02,1.09,3.38,3.42,3.49,2.91,3.85,2.79,2.56,1.92,
    2.26,2.20,2.18,2.10,2.07,2.05,1.93,1.74,1.50,1.27,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.08,0.21,0.38,0.68,0.44,0.32,0.29,0.18,0.11,0.09,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.10,0.98,2.41,3.29,3.38,4.29,3.70,3.51,2.90,
    1.94,2.00,2.08,2.16,2.15,2.14,2.12,1.95,1.48,1.39,
    0.  ,0.02,0.07,0.28,0.39,0.59,0.76,0.97,0.91,0.85,
    0.82,0.80,0.78,0.76,0.74,0.74,0.63,0.59,0.52,0.4,
    0.05,0.26,0.37,0.23,0.24,0.20,0.17,0.10,0.08,0.07,
    0.07,0.09,0.20,0.14,0.12,0.11,0.06,0.05,0.05,0.05,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.06,0.04,0.03,0.02,0.02,0.11,0.09,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.04,0.10,0.16,0.22,0.24,0.35,0.41,
    0.41,0.40,0.42,0.42,0.43,0.43,0.44,0.44,0.43,0.39,
    0.  ,0.  ,0.  ,0.07,0.31,0.48,0.64,0.66,0.83,0.98,
    0.74,0.38,0.40,0.41,0.41,0.42,0.43,0.41,0.26,0.24,
    0.  ,0.02,0.07,0.06,0.21,0.34,0.56,0.48,0.59,0.64,
    0.68,0.71,0.54,0.46,0.38,0.38,0.41,0.36,0.32,0.33,
    0.  ,0.  ,0.  ,0.10,0.26,0.38,0.67,0.41,0.51,0.48,
    0.36,0.21,0.24,0.22,0.21,0.21,0.18,0.17,0.17,0.16,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.,
    2.05,5.30,4.72,3.70,2.60,2.70,2.20,1.88,1.90,1.43,
    1.10,0.74,0.81,0.80,0.62,0.48,0.47,0.47,0.47,0.66,
    2.12,1.86,0.81,0.68,0.91,0.71,0.96,1.78,1.69,1.70,
    1.30,0.94,1.03,1.04,0.84,0.84,0.67,0.67,0.66,0.65,
    4.93,10.40,6.67,5.40,4.00,3.50,3.70,2.80,2.60,2.00,
    1.50,1.14,1.26,1.28,1.06,0.90,0.88,0.87,0.85,0.83,
    0.05,0.18,0.43,0.67,0.92,1.16,1.78,1.71,1.52,1.34,
    1.04,0.97,0.92,1.02,0.89,0.82,0.77,0.73,0.70,0.66,
    0.08,0.32,0.70,1.60,1.65,1.88,1.68,1.92,1.79,1.64,
    1.49,1.37,1.30,1.70,1.40,1.26,1.19,1.13,1.08,1.02,
    0.30,0.94,0.70,0.70,0.70,0.70,0.31,0.29,0.27,0.25,
    0.23,0.21,0.18,0.15,0.13,0.11,0.09,0.07,0.05,0.03,
    0.33,1.19,1.89,3.56,4.20,3.78,3.35,2.93,2.51,2.12,
    1.40,1.34,1.28,1.22,1.16,1.10,1.03,1.00,0.97,0.94,
    0.20,0.22,0.23,1.16,1.53,1.86,1.91,2.18,2.02,1.83,
    1.65,1.52,1.72,1.55,1.39,1.23,1.11,1.02,1.36,1.69,
    0.30,0.40,0.70,2.20,2.10,1.80,1.90,2.10,2.20,2.00,
    1.80,1.70,1.40,1.30,1.35,1.30,1.28,1.10,1.17,1.24,
    0.09,0.12,0.30,0.37,0.67,0.89,0.98,1.08,1.16,1.16,
    1.16,1.17,1.17,1.21,0.79,0.55,0.66,0.79,0.89,1.00,
    0.  ,0.  ,0.09,0.11,0.05,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.15,0.26,0.43,0.35,0.28,0.14,0.10,0.08,
    0.04,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.10,0.28,0.44,0.85,1.15,1.87,2.22,1.80,
    2.94,1.98,1.72,1.50,1.16,1.05,0.94,0.95,1.06,1.17,
    0.  ,0.02,0.04,0.12,0.60,1.60,1.30,1.50,2.10,2.40,
    2.65,2.40,2.09,1.94,1.83,1.17,1.06,0.94,0.93,0.98,
    0.  ,.003,.041,0.13,0.40,0.09,0.10,0.38,0.96,1.07,
    2.13,2.14,2.07,2.00,1.95,1.83,1.17,0.86,0.98,0.90,
    0.  ,.001,0.09,0.19,0.52,0.64,0.78,0.85,0.99,1.17,
    1.86,1.93,2.14,1.98,1.83,1.76,1.52,1.43,1.35,1.24,
    0.10,0.10,0.03,0.16,0.35,0.37,0.52,0.64,0.82,0.95,
    1.09,1.16,1.73,1.80,1.12,1.03,0.95,0.93,0.90,0.84,
    0.  ,.015,.008,0.11,0.20,0.21,0.34,0.58,0.46,0.38,
    0.65,1.17,1.58,1.61,1.73,1.82,1.74,1.63,1.23,1.00,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.
  };

  const G4int ifkn[7][7][18] = {
    1,1,7,1,1,7,7,1,1,7,7,7,1,1,7,7,7,7,3,2,
    1,1,1,3,5,1,1,7,3,5,1,1,7,7,3,5,0,0,0,1,
    2,7,3,1,2,3,3,5,1,1,3,5,3,5,0,0,0,2,2,3,
    3,1,2,7,7,3,1,2,3,5,3,7,0,0,0,0,0,0,0,2,
    2,7,3,3,1,2,7,7,7,3,0,0,0,0,0,0,0,0,0,0,
    0,0,2,2,3,3,7,7,0,3,2,1,5,1,7,2,7,0,0,0,2,2,3,3,5,3,
    1,2,7,1,2,7,7,1,2,7,7,7,1,2,7,7,7,7,1,1,
    5,1,2,3,5,1,2,3,5,7,1,2,3,5,3,5,2,2,3,1,
    1,5,7,1,1,5,3,5,1,2,3,5,7,7,0,0,0,2,2,3,
    7,1,1,5,7,7,1,1,3,5,5,7,0,0,0,0,0,0,0,2,
    2,3,7,7,1,1,7,7,7,5,0,0,0,0,0,0,0,2,2,3,
    5,3,2,2,7,7,7,3,0,1,7,2,7,2,3,1,5,0,0,0,2,2,3,5,7,3,
    3,1,7,3,1,7,7,3,1,7,7,7,3,1,7,7,7,7,3,3,
    2,3,3,2,7,3,3,5,1,7,3,5,3,5,3,1,0,0,0,3,
    5,1,3,3,5,3,3,2,3,5,7,7,3,1,0,0,0,0,0,0,
    0,3,3,2,7,7,3,5,3,3,7,2,0,0,0,0,0,0,0,0,
    0,0,0,0,3,3,2,7,7,7,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    2,2,7,2,2,7,7,2,2,7,7,7,2,2,7,7,7,7,5,2,
    1,2,2,3,5,2,2,3,5,7,2,2,3,5,7,7,0,0,0,2,
    1,7,5,2,1,5,5,3,2,2,3,5,3,5,0,0,0,1,1,5,
    5,1,2,7,7,5,1,2,3,5,5,7,0,0,0,0,0,0,0,1,
    1,7,5,5,1,2,7,7,7,5,0,0,0,0,0,0,0,0,0,0,
    0,0,1,1,5,5,7,7,0,0,0,0,0,0,0,0,0,0,0,0,1,1,5,5,3,5,
    5,1,7,5,1,7,7,5,1,7,7,7,5,1,3,5,3,5,7,7,
    2,5,1,3,5,5,1,5,3,7,5,1,7,7,7,7,3,5,2,7,
    7,7,2,7,7,7,7,2,5,1,5,3,7,7,0,0,0,3,5,2,
    7,3,5,2,7,7,2,7,7,7,7,7,0,0,0,0,0,0,0,3,
    5,3,5,2,3,5,2,7,7,7,0,0,0,0,0,0,0,0,0,0,
    0,0,3,5,3,5,2,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    2,3,7,3,2,7,7,3,2,7,7,7,3,2,3,5,3,5,7,7,
    1,3,2,3,5,3,2,3,5,7,3,2,7,7,7,7,3,5,1,7,
    7,7,1,1,7,7,7,7,3,2,3,5,7,7,0,0,0,3,5,1,
    7,3,5,1,7,7,1,7,7,7,7,7,0,0,0,0,0,0,0,1,
    3,5,3,5,3,5,1,7,7,7,0,0,0,0,0,0,0,0,0,0,
    0,0,3,5,3,5,1,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    5,2,7,5,2,7,7,5,2,7,7,7,7,7,7,7,5,2,5,5,
    1,5,1,7,5,3,5,7,5,2,3,5,3,5,5,2,0,0,0,5,
    3,2,5,3,5,5,5,1,3,5,7,7,5,2,0,0,0,0,0,0,
    0,5,5,1,7,7,3,5,5,5,1,7,0,0,0,0,0,0,0,0,
    0,0,0,0,5,5,1,7,7,7,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
  };

  vector<G4int> kinds;

  G4int l = is;

  if(l == 10) {
    l = 3;
  }
  else if(l == 4) {
    l = 1;
  }
  else if(l == 6 || l == 5) {
    l = 4;
  };

  G4int il = getIL(is, mult);
  pair<G4int, G4double> iksk = getPositionInEnergyScale1(ekin);
  G4int ik = iksk.first;
  G4double sk = iksk.second;
  G4int n;      

  switch(mult) {
  case 3:
    n = 0;
    break;
  case 4:
    n = 3;
    break;
  case 5:
    n = 7;
    break;
  case 6:
    n = 13;
    break;
  default:
    G4cout << " generateOutgoingKindsFor2toMany: mult " << mult << G4endl;
    n = 13;
  };      

  vector<G4double> sig;
  G4double stot = 0.0;

  if(l == 7 || l == 14) {
    for(G4int j = 0; j < il; j++) {
      sig.push_back(0.5 * (bsig[2][n + j][ik - 1] + bsig[3][n + j][ik - 1] + 
			   sk * (bsig[2][n + j][ik] - bsig[2][n + j][ik - 1] +
				 bsig[3][n + j][ik] - bsig[3][n + j][ik - 1])));
      stot += sig[j];
    };
  }
  else {
    for(G4int j = 0; j < il; j++) { 
      sig.push_back(bsig[l - 1][n + j][ik - 1] + 
		    sk*(bsig[l - 1][n + j][ik] - bsig[l - 1][n + j][ik - 1]));
      stot += sig[j];
    };
  };

  G4double sl = inuclRndm();

  sl *= stot;

  G4double ptot = 0.0;
  G4int ml;

  for(G4int i = 0; i < il; i++) {
    ptot += sig[i];
    if(sl <= ptot) {
      ml = i;
      break;
    };  
  }; 

  l = is;
  if(l == 14) {
    l = 5;
  }
  else if(l == 7) {
    l = 6;  
  }
  else if(is == 10) {
    l = 7;
  };

  G4int ks;      

  switch(mult) {
  case 3:
    ks = 0;
    break;
  case 4:
    ks = 3;
    break;
  case 5:
    ks = 7;
    break;
  case 6:
    ks = 12;
    break;
  default:
    G4cout << " generateOutgoingKindsFor2toMany: mult " << mult << G4endl;
    ks = 12;
  };      
  for(i = 0; i < mult; i++)  
    kinds.push_back(ifkn[l - 1][ml][i + ks]);

  return kinds;
}	

G4double G4ElementaryParticleCollider::getMomModuleFor2toMany( 
							      G4int is, 
							      G4int mult, 
							      G4int knd, 
							      G4double ekin) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::getMomModuleFor2toMany" << G4endl;
  }

  const G4double rmn[14][10][2] = {
    0.5028,  0.6305,  3.1442, -3.7333, -7.8172,  13.464,  8.1667, 
    -18.594, 1.6208,  1.9439, -4.3139, -4.6268,  12.291,  9.7879, 
    -15.288, -9.6074, 0.0,     0.0,     0.0,     0.0,     0.9348,  
    2.1801,  -10.59,  1.5163,  29.227, -16.38,  -34.55,  27.944,
    -0.2009,-0.3464,  1.3641,  1.1093, -3.403,  -1.9313,  3.8559,
    1.7064,  0.0,     0.0,     0.0,     0.0,    -0.0967, -1.2886, 
    4.7335, -2.457, -14.298,   15.129,  17.685, -23.295,  0.0126,  
    0.0271, -0.0835, -0.1164,  0.186,   0.2697, -0.2004, -0.3185, 
    0.0,     0.0,     0.0,     0.0,    -0.025,   0.2091, -0.6248,  
    0.5228, 2.0282, -2.8687,  -2.5895,  4.2688, -0.0002, -0.0007,  
    0.0014,  0.0051, -0.0024,-0.015,    0.0022,  0.0196,  0.0,     
    0.0,     0.0,     0.0,     1.1965,  0.9336, -0.8289, -1.8181,
    1.0426,  5.5157, -1.909,  -8.5216,  1.2419,  1.8693, -4.3633,
    -5.5678, 13.743,  14.795, -18.592, -16.903,  0.0,     0.0,
    0.0,     0.0,     0.287,   1.7811, -4.9065, -8.2927,  16.264,
    20.607,-19.904, -20.827,  -0.244,  -0.4996,  1.3158,  1.7874,
    -3.5691,-4.133,   4.3867,  3.8393,  0.0,     0.0,     0.0,
    0.0,    -0.2449, -1.5264,  2.9191,  6.8433, -9.5776, -16.067,
    11.938, 16.845,   0.0157,  0.0462, -0.0826, -0.1854,  0.2143,
    0.4531, -0.2585, -0.4627,  0.0,     0.0,     0.0,     0.0,
    0.0373,  0.2713, -0.422,  -1.1944,  1.3883,  2.7495, -1.7476,
    -2.9045,-0.0003, -0.0013,  0.0014,  0.0058, -0.0034, -0.0146,
    0.0039,  0.0156,  0.0,     0.0,     0.0,     0.0,     0.0,    
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,    
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     
    0.0,     0.1451,  0.0929,  0.1538,  0.1303,  0.0,     0.0,     
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     
    0.4652,  0.5389,  0.2744,  0.4071,  0.0,     0.0,     0.0,
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0, 
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,    -0.033,
    -0.0545,-0.0146, -0.0288,  0.0,     0.0,     0.0,     0.0,
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
    0.0,     0.0,     0.0,     0.0,     0.0,     0.6296,  0.1491,
    0.8381,  0.1802,  0.0,     0.0,     0.0,     0.0,     0.0,
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
    0.0,     0.0,     0.0,     0.0,     0.1787,  0.385,   0.0086,
    0.3302,  0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
    0.0,     0.0,     0.0,     0.0,     0.0,     0.0,     0.0,
    0.0,     0.0,     0.0,    -0.0026, -0.0128,  0.0033, -0.0094
  };

  G4double S = inuclRndm();
  G4double PS = 0.0;
  G4double PR = 0.0;
  G4double PQ = 0.0;
  G4int KM = 2;
  G4int IL = 4;
  G4int JK = 4;
  G4int JM = 2;
  G4int IM = 3;

  if(is == 1 || is == 2 || is == 4) KM = 1;
  if(mult == 3) { IM = 0; IL = 0; };
  if(knd == 1 || knd == 2) JK = 0;
  for(G4int i = 0; i < 4; i++) {

    G4double V = 0.0;

    for(G4int k = 0; k < 4; k++) V += rmn[k + JK][i + IL][KM - 1] * pow(ekin, k);
    PR += V * pow(S, i);
    PQ += V;
  };  

  if(knd == 1 || knd == 2) JM = 1;
  for(G4int m = 0; m < 3; m++) PS += rmn[8 + IM + m][7 + JM][KM - 1] * pow(ekin, m);

  G4double PRA = PS * sqrt(S) * (PR + (1 - PQ) * pow(S, 4));

  return fabs(PRA);
}

vector<G4double> G4ElementaryParticleCollider::
particleSCMmomentumFor2to3(
			   G4int is, 
			   G4int knd, 
			   G4double ekin, 
			   G4double pmod) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::particleSCMmomentumFor2to3" << G4endl;
  }

  const G4double abn[4][4][4] = {
    {{0.0856,  0.0716,  0.1729,  0.0376},  {5.0390,  3.0960,  7.1080,  1.4331},
    {-13.782, -11.125, -17.961, -3.1350},  {14.661,  18.130,  16.403,  6.4864}},
    {{0.0543,  0.0926, -0.1450,  0.2383}, {-9.2324, -3.2186, -13.032,  1.8253},
    {36.397,  20.273,  41.781,  1.7648}, {-42.962, -33.245, -40.799, -16.735}},
    {{-0.0511, -0.0515,  0.0454, -0.1541}, {4.6003,  0.8989,  8.3515, -1.5201},
    {-20.534, -7.5084, -30.260, -1.5692},  {27.731,  13.188,  32.882,  17.185}},
    {{0.0075,  0.0058, -0.0048,  0.0250}, {-0.6253, -0.0017, -1.4095,  0.3059},
     {2.9159,  0.7022,  5.3505,  0.3252}, {-4.1101, -1.4854, -6.0946, -3.5277}} 
  };


  /*
  const G4double abn[4][4][4] = {
    0.0856,  0.0716,  0.1729,  0.0376,  5.0390,  3.0960,  7.1080,  1.4331,
    -13.782, -11.125, -17.961, -3.1350,  14.661,  18.130,  16.403,  6.4864,
    0.0543,  0.0926, -0.1450,  0.2383, -9.2324, -3.2186, -13.032,  1.8253,
    36.397,  20.273,  41.781,  1.7648, -42.962, -33.245, -40.799, -16.735,
    -0.0511, -0.0515,  0.0454, -0.1541,  4.6003,  0.8989,  8.3515, -1.5201,
    -20.534, -7.5084, -30.260, -1.5692,  27.731,  13.188,  32.882,  17.185,
    0.0075,  0.0058, -0.0048,  0.0250, -0.6253, -0.0017, -1.4095,  0.3059,
    2.9159,  0.7022,  5.3505,  0.3252, -4.1101, -1.4854, -6.0946, -3.5277 
  };
  */

  const G4int itry_max = 100;
  G4double ct = 2.0;
  G4int K = 3;
  G4int J = 1;

  if(is == 1 || is == 2 || is == 4) K = 1;
  if(knd == 1 || knd == 2) J = 0;

  G4int itry = 0;

  while(fabs(ct) > 1.0 && itry < itry_max) {
    itry++;

    G4double S = inuclRndm();
    G4double U = 0.0;
    G4double W = 0.0;

    for(G4int l = 0; l < 4; l++) {

      G4double V = 0.0;

      for(G4int m = 0; m < 4; m++) {
	V += abn[m][l][K + J - 1] * pow(ekin, m);
      };
      U += V;
      W += V * pow(S, l);
    };  
    ct = 2.0 * sqrt(S) * (W + (1.0 - U) * pow(S, 4)) - 1.0;
  };
  if(itry == itry_max) {

    if(verboseLevel > 2){
      G4cout << " particleSCMmomentumFor2to3 -> itry = itry_max " << itry << G4endl;
    }

    ct = 2.0 * inuclRndm() - 1.0;
  };

  G4double pt = pmod * sqrt(1.0 - ct * ct);
  G4double phi = randomPHI();

  vector<G4double> mom(4);

  mom[1] = pt * cos(phi);
  mom[2] = pt * sin(phi);
  mom[3] = pmod * ct;
  
  return mom;  
}
	      
G4bool G4ElementaryParticleCollider::reChargering(G4double ekin, 
						  G4int is) const {
  
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::reChargering" << G4endl;
  }

  const G4double ali[31] = {
    3.9, 6.6, 13.2, 42.6, 36.6, 26.2, 17.9, 13.5, 11.3, 10.7,
    9.5, 6.4,  1.9, 0.8,  0.56, 0.26, 0.19, 0.19, 0.18, 0.1,
    0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 ,0.0
  };

  const G4double asig[4][6][31] = {
    1.00,1.00,1.00,1.00,1.00,1.00,24.3,24.1,24.0,26.3,
    28.6,24.8,19.9,19.2,17.4,15.3,13.5,12.3,11.9,10.4,
    11.8,11.4,11.0,10.8,10.9,11.7,11.4,10.2,11.0,11.0,
    9.0, 0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,1.45,2.90,
    4.10,5.30,22.0,21.2,19.9,14.9,12.6,12.3,11.3,10.8,
    9.50,8.27,7.20,6.70,6.25,6.04,5.89,5.70,5.60,5.05,
    4.17,4.00,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,2.27,4.28,7.81,8.03,8.11,7.90,7.82,
    7.61,7.47,7.07,7.66,7.05,6.71,6.38,6.36,6.37,6.57,
    6.01,5.48,6.80,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,1.20,2.85,3.70,4.81,5.33,
    7.74,6.91,6.94,7.57,7.21,7.11,7.10,6.93,6.79,6.71,
    6.55,6.55,6.15,8.50,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,.005,0.54,0.74,
    0.86,0.91,1.10,1.16,1.36,1.40,1.43,1.47,1.47,1.43,
    1.38,1.38,1.63,1.36,2.80,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,34.0,46.2,46.9,45.2,47.1,
    42.3,41.8,41.2,41.6,41.6,41.0,43.0,42.4,40.0,39.9,
    39.8,42.0,40.0,39.8,39.6,38.7,
    1.00,1.00,1.00,1.00,1.00,1.00,33.0,31.3,29.5,27.8,
    14.6,16.0,17.5,18.3,19.4,18.7,15.6,14.8,13.6,12.5,
    12.2,11.9,11.4,11.2,10.1,9.62,8.41,7.14,7.09,5.04,
    10.2,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,2.00,4.00,
    6.50,19.5,20.8,19.8,18.6,17.7,14.4,13.5,10.4,10.1,
    12.0,8.87,8.51,8.49,9.20,8.29,7.43,8.20,4.69,4.57,
    4.06,4.10,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.68,2.76,3.85,8.35,12.7,12.3,10.8,12.0,
    10.9,10.2,10.6,12.2,11.8,11.0,10.4,10.5,11.0,10.6,
    11.8,12.5,13.2,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.22,0.52,0.64,1.19,1.52,1.75,1.51,
    2.04,1.85,1.70,1.92,1.66,1.74,1.50,1.39,1.35,1.41,
    1.48,1.43,1.35,3.20,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.14,0.24,0.30,0.46,
    0.85,1.40,1.54,1.52,1.47,1.48,1.49,1.42,1.39,1.37,
    1.22,1.19,0.93,0.  ,2.10,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,35.0,40.0,42.4,42.3,41.0,
    40.9,40.4,39.8,35.0,33.6,41.2,41.0,41.1,41.2,41.2,
    39.6,36.0,36.0,36.2,0.  ,40.2,
    1.00,1.00,1.00,1.00,1.00,1.00,60.0,38.0,30.6,24.0,
    18.5,12.8,13.6,9.15,8.20,7.80,7.10,6.40,5.81,5.85,
    5.50,5.33,5.40,5.50,4.90,5.02,5.00,4.98,4.96,4.96,
    4.50,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.38,0.75,
    1.28,2.26,11.8,11.1,7.56,6.04,5.68,4.15,3.21,2.12,
    1.66,1.54,1.53,1.47,1.33,1.39,1.35,1.29,1.23,1.13,
    1.06,0.70,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.05,2.19,6.02,7.36,6.97,5.83,6.53,5.09,
    4.24,3.24,3.31,3.11,2.91,2.78,2.72,2.68,2.48,2.27,
    2.02,1.77,4.40,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.06,0.35,0.55,0.77,0.89,0.80,0.86,
    0.97,0.90,0.86,0.84,0.84,0.83,0.81,0.79,0.82,0.76,
    0.88,0.89,0.89,3.00,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,0.02,0.07,0.33,0.92,1.39,
    2.11,1.81,2.39,2.60,2.19,1.70,1.60,0.68,1.43,1.46,
    1.46,1.37,1.16,1.09,2.60,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,18.9,27.2,34.9,29.1,30.8,
    29.6,28.2,27.5,26.9,26.3,25.9,25.6,25.2,26.1,25.5,
    25.4,25.3,25.1,24.9,24.8,24.1,
    5.90,9.40,24.5,62.6,65.3,41.3,29.3,24.3,22.7,22.9,
    23.2,28.4,11.7,10.1,8.30,7.16,6.49,6.36,6.60,5.84,
    5.30,4.50,3.90,4.40,4.74,.794,.824,.714,0.59,0.  ,
    4.60,0.  ,0.  ,0.  ,0.  ,0.10,0.40,2.70,3.50,5.30,
    6.60,9.10,17.6,12.2,9.78,7.51,6.91,6.86,6.46,6.19,
    5.13,3.90,2.82,3.10,3.12,2.52,2.22,2.02,2.01,1.98,
    2.14,1.20,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.76,2.63,3.72,6.53,7.47,7.94,7.12,6.85,
    6.09,5.35,4.12,3.85,3.68,4.09,3.58,3.29,3.08,2.93,
    2.80,2.65,3.30,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.59,0.74,1.47,4.10,4.78,4.90,5.07,
    5.50,5.48,5.03,4.65,4.39,4.06,3.53,3.08,3.05,2.91,
    3.42,3.93,3.93,4.10,0.  ,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.01,.007,0.03,.099,.251,.376,
    .419,.582,.755,.777,1.13,1.08,1.13,1.08,.962,.866,
    .738,.674,.645,.613,1.30,0.  ,0.  ,0.  ,0.  ,0.  ,
    0.  ,0.  ,0.  ,0.  ,0.  ,31.3,46.0,30.0,35.7,33.4,
    31.6,30.4,29.6,28.9,28.5,28.1,27.5,31.0,27.7,27.8,
    26.1,25.2,6.92,6.70,0.  ,25.7
  };

  G4bool rech = false;

  if(is == 6 || is == 5 || is == 7 || is == 14) {
    pair<G4int, G4double> iksk = getPositionInEnergyScale2(ekin);
    G4int ik = iksk.first;
    G4double sk = iksk.second;
    G4double chrg;

    if(ik == 30 && sk == 1.0) {
      chrg = 1.0;    
    }
    else {
      chrg = 1.0 - (ali[ik - 1] + sk * (ali[ik] - ali[ik - 1])) /
	(asig[3][0][ik - 1] + sk * (asig[3][0][ik] - asig[3][0][ik - 1]));
    }; 
    rech = inuclRndm() > chrg; 
  };

  return rech;
}

pair<G4double, G4double> G4ElementaryParticleCollider::
adjustIntervalForElastic(
			 G4double ekin, 
			 G4double ak, 
			 G4double ae, 
			 G4int k, 
			 G4int l, 
			 const vector<G4double>& ssv,
			 G4double st) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::adjustIntervalForElastic" << G4endl;
  }

  const G4double ang[4][4][13] = {
    2.7404,-30.853, 0.1026,-0.3829, 0.2499, 3.9025, 19.402, 0.1579,
    0.3153,-17.953, 0.4217, 0.1499, 0.5369,-9.6998, 106.24,-1.0542,
    3.7587, 32.028,-91.126,-224.46, 2.9671,-7.4981, 109.72, 147.05,
    2.8753,-13.216, 10.400,-129.39, 11.389,-6.5144,-118.82, 323.73,
    747.33,-5.5251, 43.295,-239.54,-653.35,-5.3078, 81.011, 2.3882,
    54.339,-16.638, 6.7740, 150.99,-400.48,-935.70, 6.8925,-76.460,
    228.26, 915.07, 6.2233,-142.85,-7.5137, 19.465,-0.4961, 103.81,
    -2.6994,-20.619,-44.180,-7.0218,-6.5373, 91.968,-3.5198,-5.9558,
    -10.550, 44.096,-68.102, 11.800,-272.82,-460.45, 491.70, 471.94,
    -205.34, 193.07,-519.63,-260.19,-162.03, 296.29,-74.379, 96.358,
    -90.857, 477.59, 1895.9,-1715.5,-1485.6, 569.51,-1018.1, 1126.6,
    1225.0, 430.79,-1695.7, 46.038,-56.827, 164.76,-512.22,-2519.0,
    2114.3, 1805.5,-898.58, 1742.6,-1074.0,-1748.1,-625.48, 2893.5,
    7.5479,-3.4831, 1.5437,-1788.2, 16.268, 33.004, 31.567, 134.96,
    46.864,-132.70, 3.6373, 128.75, 69.621,-39.274, 12.341,-33.769,
    4305.2, 2138.4,-766.84,-301.76, 4872.2,-1303.0, 741.12, 155.92,
    3140.2,-1924.5, 64.835,-18.592, 251.92,-7931.4,-9126.2, 2700.3,
    907.63,-14674., 6729.1,-1600.0,-752.01,-7918.9, 10620., 41.609,
    12.024,-450.71, 9347.1, 12431.,-3352.5,-1077.3, 23924.,-11075.,
    1524.9, 1079.6, 10983.,-17468.,-1.8369, 0.1894,-1.2021, 7147.5,
    -29.654,-16.367,-6.8648,-821.16,-95.192, 58.598,-0.7804,-851.61,
    -138.65, 8.6911,-0.6788, 0.2534,-3339.5,-3182.3, 373.94, 60.476,
    -32586., 2637.3,-318.74,-30.563,-18780., 3928.1,-13.060, 1.0665,
    -186.58,-4139.2, 13944.,-1320.2,-175.20,100980.,-12857., 677.51,
    147.95, 44607.,-20293., 7.1880,-0.7291, 332.54,-4436.4,-19342.,
    1642.3, 203.81,-165530.,20294.,-640.11,-212.50,-58790., 32058.
  };
  
  const G4int itry_max = 100;
  const G4double small = 1.0e-4;

  G4double a = 1.0;
  G4double b = 0.0;
  G4int adj_type = 0;
  G4double s1;
  G4double s2;

  if(k == 1) {
    adj_type = 1;
    s1 = 0.0;
    s2 = 0.5;
  }
  else if(k == 2) {
    if(l != 2) {
      adj_type == 0;
    }
    else { 
      adj_type = 1;
      s1 = 0.0;
      s2 = 0.67;
    };        
  }
  else {
    if(ekin < 0.32) {
      adj_type == 0;
    }
    else {
      adj_type = 2;
      s1 = 0.1813;
      s2 = 1.0;
    };     
  };

  if(adj_type > 0) {

    G4int itry = 0;
    G4double su;
    G4double ct;
  
    if(adj_type == 1) {

      G4double s2_old = s2;
      G4double s1c = s1;
      G4double s2_new;

      while(itry < itry_max) {
	itry++;
	s2_new = 0.5 * (s2_old + s1c);
	su = 0.0;      
	for(G4int i = 0; i < 4; i++) su += ssv[i] * pow(s2_new, i);
	ct = ak * sqrt(s2_new) * (su + (1.0 - st) * pow(s2_new, 4)) + ae;
	if(ct > 1.0) {
	  s2_old = s2_new;
	}
	else {
	  if(1.0 - ct < small) {
	    break;
	  }
	  else {
	    s1c = s2_new;
	  };   
	}; 
      };
      a = s2_new - s1;
      b = s1;
    }
    else {

      G4double s1_old = s1;
      G4double s2c = s2;
      G4double s1_new;

      while(itry < itry_max) {
	itry++;
	s1_new = 0.5 * (s1_old + s2c);
	su = 0.0;      
	for(G4int i = 0; i < 4; i++) su += ssv[i] * pow(s1_new, i);
	ct = ak * sqrt(s1_new) * (su + (1.0 - st) * pow(s1_new, 4)) + ae;
	if(ct < -1.0) {
	  s1_old = s1_new;
	}
	else {
	  if(1.0 + ct < small) {

	    break;
	  }
	  else {
	    s2c = s1_new;
	  };   
	}; 
      };
      a = s2 - s1_new;
      b = s1_new;
    }; 
    if(itry == itry_max) {

      if(verboseLevel > 2){
	G4cout << " in adjustIntervalForElastic: " << itry_max << G4endl;
	G4cout << " e " << ekin << " ak " << ak << " ae " << ae << G4endl << " k " << k
	       << " is " << l << " adj_type " << adj_type << G4endl; 
      }

      a = 1.0;
      b = 0.0;
    };
  };

  return pair<G4double, G4double>(a, b);
}  

vector<G4double> G4ElementaryParticleCollider::
particleSCMmomentumFor2to2(
			   G4int is, 
			   G4int kw, 
			   G4double ekin, 
			   G4double pscm) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::particleSCMmomentumFor2to2" << G4endl;
  }

  const G4double ang[4][4][13] = {
    2.7404,-30.853, 0.1026,-0.3829, 0.2499, 3.9025, 19.402, 0.1579,
    0.3153,-17.953, 0.4217, 0.1499, 0.5369,-9.6998, 106.24,-1.0542,
    3.7587, 32.028,-91.126,-224.46, 2.9671,-7.4981, 109.72, 147.05,
    2.8753,-13.216, 10.400,-129.39, 11.389,-6.5144,-118.82, 323.73,
    747.33,-5.5251, 43.295,-239.54,-653.35,-5.3078, 81.011, 2.3882,
    54.339,-16.638, 6.7740, 150.99,-400.48,-935.70, 6.8925,-76.460,
    228.26, 915.07, 6.2233,-142.85,-7.5137, 19.465,-0.4961, 103.81,
    -2.6994,-20.619,-44.180,-7.0218,-6.5373, 91.968,-3.5198,-5.9558,
    -10.550, 44.096,-68.102, 11.800,-272.82,-460.45, 491.70, 471.94,
    -205.34, 193.07,-519.63,-260.19,-162.03, 296.29,-74.379, 96.358,
    -90.857, 477.59, 1895.9,-1715.5,-1485.6, 569.51,-1018.1, 1126.6,
    1225.0, 430.79,-1695.7, 46.038,-56.827, 164.76,-512.22,-2519.0,
    2114.3, 1805.5,-898.58, 1742.6,-1074.0,-1748.1,-625.48, 2893.5,
    7.5479,-3.4831, 1.5437,-1788.2, 16.268, 33.004, 31.567, 134.96,
    46.864,-132.70, 3.6373, 128.75, 69.621,-39.274, 12.341,-33.769,
    4305.2, 2138.4,-766.84,-301.76, 4872.2,-1303.0, 741.12, 155.92,
    3140.2,-1924.5, 64.835,-18.592, 251.92,-7931.4,-9126.2, 2700.3,
    907.63,-14674., 6729.1,-1600.0,-752.01,-7918.9, 10620., 41.609,
    12.024,-450.71, 9347.1, 12431.,-3352.5,-1077.3, 23924.,-11075.,
    1524.9, 1079.6, 10983.,-17468.,-1.8369, 0.1894,-1.2021, 7147.5,
    -29.654,-16.367,-6.8648,-821.16,-95.192, 58.598,-0.7804,-851.61,
    -138.65, 8.6911,-0.6788, 0.2534,-3339.5,-3182.3, 373.94, 60.476,
    -32586., 2637.3,-318.74,-30.563,-18780., 3928.1,-13.060, 1.0665,
    -186.58,-4139.2, 13944.,-1320.2,-175.20,100980.,-12857., 677.51,
    147.95, 44607.,-20293., 7.1880,-0.7291, 332.54,-4436.4,-19342.,
    1642.3, 203.81,-165530.,20294.,-640.11,-212.50,-58790., 32058.
  };

  const G4int itry_max = 100;
  const G4double ct_cut = 0.9999;
  const G4double huge = 60.0;
  G4int k = getElasticCase(is, kw, ekin);
  G4double ae = -1.0;
  G4double ak = 2.0;

  if(k == 1) { 
    if(is != 2) { ak = 1.0; ae = 0.0;};
  }
  else if(k == 2) {
    if(is != 2) { ak = 0.5; ae = 0.0; };
  };    

  G4double ct = 2.0;
  G4double ab;
  G4double ac;
  G4double ad;
  G4int itry = 0;
  
  if(k == 14) {
    ab = 7.5;
    if(is == 1 || is == 2 || is == 4) ab = 8.7;
    ac = -2.0 * ab * pscm * pscm;
    ad = 2.0 * ac;
    if(ad < -huge) {
      ad = exp(ad);
    }
    else {
      ad = exp(-huge);
    };   
    while(fabs(ct) > ct_cut && itry < itry_max) {
      itry++;
      ct = 1.0 - log(inuclRndm() * (1.0 - ad) + ad) / ac;
    };      
  }
  else if(k == 0) {
    ct = 2.0 * inuclRndm() - 1;
  }
  else { 
    G4int k1 = k - 1;
    // first set all coefficients
    vector<G4double> ssv(4);
    G4double st = 0.0;

    for(G4int i = 0; i < 4; i++) {

      G4double ss = 0.0;

      for(G4int m = 0; m < 4; m++) ss += ang[m][i][k1] * pow(ekin, m);
      st += ss;
      ssv[i] = ss;
    };

    G4double a = 1.0;
    G4double b = 0.0;

    if(k <= 3) {
      pair<G4double, G4double> ab = adjustIntervalForElastic(ekin, ak, ae, k, is, ssv, st);

      a = ab.first;
      b = ab.second;
    };
    while(fabs(ct) > ct_cut && itry < itry_max) {
      itry++;
      G4double mrand = a * inuclRndm() + b;

      G4double su = 0.0;

      for(G4int i = 0; i < 4; i++) su += ssv[i] * pow(mrand, i);
      ct = ak * sqrt(mrand) * (su + (1.0 - st) * pow(mrand, 4)) + ae;
    }; 
  };

  if(itry == itry_max) {

    if(verboseLevel > 2){
      G4cout << " particleSCMmomentumFor2to2 -> itry = itry_max " << itry << G4endl;
    }

    ct = 2.0 * inuclRndm() - 1.0;
  };

  G4double pt = pscm * sqrt(1.0 - ct * ct);
  G4double phi = randomPHI();
  vector<G4double> mom(4);

  mom[1] = pt * cos(phi);
  mom[2] = pt * sin(phi);
  mom[3] = pscm * ct;
  
  return mom;  
}

G4int G4ElementaryParticleCollider::getElasticCase(G4int is, 
						   G4int kw, 
						   G4double ekin) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::getElasticCase" << G4endl;
  }

  G4int l = is;
  G4int k = 0; // isotropic

  if(l == 4) {
    l = 1;
  } 
  else if(l == 10 || l == 7 || l == 14) {
    l = 3;
  }
  else if(l == 5 || l == 6) {
    l = 4;
  };
  if(l < 3) { // nucleon nucleon
    if(ekin > 2.8) {
      k = 2;
      if(ekin > 10.0) k = 14;
    }
    else {    
      if(l == 1) { // PP or NN
        if(ekin > 0.46) k = 1;
      }
      else {
        k = 3;
	if(ekin >= 0.97) k = 1;
      }; 
    };  	 
  }
  else { // pi nucleon
    if(l == 3) { // pi+ P, pi- N, pi0 P, pi0 N
      k = 8;
      if(ekin > 0.08) k = 9;
      if(ekin > 0.3) k = 10;
      if(ekin > 1.0) k = 11;
      if(ekin > 2.4) k = 14;
    }
    else { // pi- P, pi+ N
      if(kw == 1) {
        k = 4;
        if(ekin > 0.08) k = 5;
        if(ekin > 0.3) k = 6;
        if(ekin > 1.0) k = 7;
        if(ekin > 2.4) k = 14;
      }
      else {
        k = 12;
	if(ekin > 0.08) k = 13;
        if(ekin > 0.3) k = 6;
        if(ekin > 1.0) k = 7;
        if(ekin > 2.4) k = 14;
      }; 
    };   
  };     

  return k;
}

vector<G4InuclElementaryParticle> G4ElementaryParticleCollider:: 
generateSCMpionAbsorption(G4double etot_scm,
			  G4InuclElementaryParticle* particle1,
			  G4InuclElementaryParticle* particle2) const { 

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::generateSCMpionAbsorption" << G4endl;
  }

  // generate nucleons momenta for pion absorption
  // the nucleon distribution assumed to be isotropic in SCM

  G4InuclElementaryParticle dummy;
  vector<G4InuclElementaryParticle> particles;
  vector<G4int> particle_kinds;
  G4int type1 = particle1->type();
  G4int type2 = particle2->type();

  // generate kinds
  if(type1 == 3) {
    if(type2 == 111) { // pi+ + PP -> ? 

      G4cout << " pion absorption: pi+ + PP -> ? " << G4endl;

      return particles;
    }
    else if(type2 == 112) { // pi+ + PN -> PP
      particle_kinds.push_back(1);
      particle_kinds.push_back(1);
    }
    else if(type2 == 122) { // pi+ + NN -> PN
      particle_kinds.push_back(1);
      particle_kinds.push_back(2);
    };     
  }
  else if(type1 == 5) { 
    if(type2 == 111) { // pi- + PP -> PN 
      particle_kinds.push_back(1);
      particle_kinds.push_back(2);
    }
    else if(type2 == 112) { // pi- + PN -> NN
      particle_kinds.push_back(2);
      particle_kinds.push_back(2);
    }
    else if(type2 == 122) { // pi- + NN -> ?

      G4cout << " pion absorption: pi- + NN -> ? " << G4endl;

      return particles;
    };     
  }
  else if(type1 == 7) {
    if(type2 == 111) { // pi0 + PP -> PP 
      particle_kinds.push_back(1);
      particle_kinds.push_back(1);
    }
    else if(type2 == 112) { // pi0 + PN -> PN
      particle_kinds.push_back(1);
      particle_kinds.push_back(2);
    }
    else if(type2 == 122) { // pi0 + NN -> ?
      particle_kinds.push_back(2);
      particle_kinds.push_back(2);
    };     
  };
    
  G4double m1 = dummy.getParticleMass(particle_kinds[0]);

  m1 *= m1;

  G4double m2 = dummy.getParticleMass(particle_kinds[1]);

  m2 *= m2;	 

  G4double a = 0.5 * (etot_scm * etot_scm - m1 - m2);

  G4double pmod = sqrt((a * a - m1 * m2) / (m1 + m2 + 2.0 * a));
  vector<G4double> mom(4);
  pair<G4double, G4double> COS_SIN = randomCOS_SIN();
  G4double FI = randomPHI();
  G4double pt = pmod * COS_SIN.second;

  mom[1] = pt * cos(FI);
  mom[2] = pt * sin(FI);
  mom[3] = pmod * COS_SIN.first;

  vector<G4double> mom1 = mom;

  for(G4int i = 1; i < 4; i++) mom1[i] *= -1.0;
  particles.push_back(G4InuclElementaryParticle(mom , particle_kinds[0]));
  particles.push_back(G4InuclElementaryParticle(mom1, particle_kinds[1]));

  return particles;
}
