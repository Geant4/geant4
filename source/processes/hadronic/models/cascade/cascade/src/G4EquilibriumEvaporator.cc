#define RUN

#include "G4EquilibriumEvaporator.hh"
#include "G4InuclNuclei.hh"
#include "G4LorentzConvertor.hh"

G4EquilibriumEvaporator::G4EquilibriumEvaporator()
  : verboseLevel(1) {
  
  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::G4EquilibriumEvaporator" << G4endl;
  }
}

G4CollisionOutput G4EquilibriumEvaporator::collide(G4InuclParticle* bullet,
						   G4InuclParticle* target) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::collide" << G4endl;
  }

  // simple implementation of the equilibium evaporation a la Dostrowski

  const G4double huge = 50.0;
  const G4double small = -50.0;
  const G4double one_third = 1.0 / 3.0;
  const G4double two_thirds = 2.0 / 3.0;
  const G4double prob_cut_off = 1.0e-15;
  const G4double Q1[6] = { 0.0, 0.0, 2.23, 8.49, 7.72, 28.3 };
  const G4double AN[6] = { 1.0, 1.0, 2.0, 3.0, 3.0, 4.0 };
  const G4double Q[6] = { 0.0, 1.0, 1.0, 1.0, 2.0, 2.0 };
  const G4double G[6] = { 2.0, 2.0, 6.0, 6.0, 6.0, 4.0 };
  const G4double BE = 0.0063;
  const G4double fisssion_cut = 1000.0;

#ifdef RUN
  const G4double cut_off_energy = 0.1;
#else
  const G4double cut_off_energy = 100000.0;
#endif

  const G4double BF = 0.0242;
  const G4int itry_max = 1000;
  const G4int itry_global_max = 1000;
  const G4double small_ekin = 1.0e-6;
  const G4int itry_gam_max = 100;

  G4std::vector<G4double> W(8);
  G4std::vector<G4double> A1(6);
  G4std::vector<G4double> Z1(6);
  G4std::vector<G4double> u(6);
  G4std::vector<G4double> V(6);
  G4std::vector<G4double> TM(6);

  G4double coul_coeff;
  G4CollisionOutput output;

  if(G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target)) {

    G4double A = nuclei_target->getA();
    G4double Z = nuclei_target->getZ();
    G4std::vector<G4double> PEX = nuclei_target->getMomentum();
    G4double EEXS = nuclei_target->getExitationEnergy();

    if(EEXS < 0.0) cout << " after noeq: eexs " << EEXS << endl;

    G4InuclElementaryParticle dummy(small_ekin, 1);
    G4LorentzConvertor toTheNucleiSystemRestFrame;

    toTheNucleiSystemRestFrame.setBullet(dummy.getMomentum(), dummy.getMass());

    G4std::vector<G4double> ppout(4, 0.0);
  
    if(timeToBigBang(A, Z, EEXS)) {

      G4cout << " big bang in eql start " << G4endl;

      return theBigBanger->collide(0, target);
    }
    else {     
      if(A >= 100.0) {
	coul_coeff = 1.4;
      }
      else {
	coul_coeff = 1.2;
      };   
   
      G4InuclNuclei dummy_nuc;
      G4double EEXS_new;
      G4std::vector<G4double> pin = PEX;

      pin[0] += 0.001 * EEXS;

      G4bool try_again = true;  
      G4bool fission_open = true;
      G4double nuc_mass;  
      G4int itry_global = 0;

      while(try_again && itry_global < itry_global_max) {
	itry_global++;

	if(verboseLevel > 2){
	  G4cout << " A " << A << " Z " << Z << " EEXS " << EEXS << G4endl;
	}

	nuc_mass = dummy_nuc.getNucleiMass(A, Z); 
	PEX[0] = sqrt(PEX[1] * PEX[1] + 
		      PEX[2] * PEX[2] + 
		      PEX[3] * PEX[3] +
		      nuc_mass * nuc_mass);  	
	toTheNucleiSystemRestFrame.setTarget(PEX, nuc_mass);
	toTheNucleiSystemRestFrame.toTheTargetRestFrame();

	if(timeToBigBang(A, Z, EEXS)) { // big bang
      
	  if(verboseLevel > 2){
	    G4cout << " big bang in eql step " << G4endl;
	  }

	  G4InuclNuclei nuclei(PEX, A, Z);        

	  nuclei.setExitationEnergy(EEXS);      

	  G4CollisionOutput explosion = theBigBanger->collide(0, &nuclei);

	  output.addOutgoingParticles(explosion.getOutgoingParticles());

	  return output;	
	}
	else { // normal chain
	  if(EEXS > cut_off_energy) { 

	    G4double E0 = getE0(A);
	    G4double parlev = getPARLEVDEN(A, Z);
	    G4double u1 = parlev * A;
	    G4std::pair<G4std::vector<G4double>, G4std::vector<G4double> > parms = paraMaker(Z);
	    G4std::vector<G4double> AK = parms.first;
	    G4std::vector<G4double> CPA = parms.second;
	    G4double DM0 = bindingEnergy(A, Z);   

	    G4int i(0);
	    for(i = 0; i < 6; i++) {
	      A1[i] = A - AN[i];
	      Z1[i] = Z - Q[i];
	      u[i] = parlev * A1[i];
	      TM[i] = -0.1;
	      if(goodRemnant(A1[i], Z1[i])) {

		G4double QB = DM0 - bindingEnergy(A1[i], Z1[i]) - Q1[i];

		V[i] = coul_coeff * Z * Q[i] * AK[i] / (1.0 + EEXS / E0) /
		  (pow(A1[i], one_third) + pow(AN[i], one_third));
		TM[i] = EEXS - QB - V[i] * A / A1[i];  
	      };
	    }; 
      
	    G4double ue = 2.0 * sqrt(u1 * EEXS);
	    G4double prob_sum = 0.0;
	
	    if(TM[0] > cut_off_energy) {

	      G4double AL = getAL(A);

	      W[0] = BE * pow(A1[0], two_thirds) * G[0] * AL;

	      G4double TM1 = 2.0 * sqrt(u[0] * TM[0]) - ue;

	      if(TM1 > huge) {
		TM1 = huge;
	      }
	      else if(TM1 < small) {
		TM1 = small;
	      };
	      W[0] = W[0] * exp(TM1); 	     
	      prob_sum += W[0];
	    }
	    else {
	      W[0] = 0.0;
	    }; 
      
	    for(i = 1; i < 6; i++) {
	      if(TM[i] > cut_off_energy) {
		W[i] = BE * pow(A1[i], two_thirds) * G[i] * (1.0 + CPA[i]);

		G4double TM1 = 2.0 * sqrt(u[i] * TM[i]) - ue;

		if(TM1 > huge) {
		  TM1 = huge;
		}
		else if(TM1 < small) {
		  TM1 = small;
		};
		W[i] = W[i] * exp(TM1); 	     
		prob_sum += W[i];
	      }
	      else {
		W[i] = 0.0;
	      }; 
	    };
	    //         fisson part
	    W[6] = 0.0;
	    if(A >= 100.0 && fission_open) {

	      G4double X2 = Z * Z / A;
	      G4double X1 = 1.0 - 2.0 * Z / A; 
	      G4double X = 0.019316 * X2 / (1.0 - 1.79 * X1 * X1);
	      G4double EF = EEXS - getQF(X, X2, A, Z, EEXS);
	  
	      if(EF > 0.0) {

		G4double AF = u1 * getAF(X, A, Z, EEXS);
		G4double TM1 = 2.0 * sqrt(AF * EF) - ue;

		if(TM1 > huge) {
		  TM1 = huge;
		}
		else if(TM1 < small) {
		  TM1 = small;
		};
		W[6] = BF * exp(TM1);
		if(W[6] > fisssion_cut * W[0]) W[6] = fisssion_cut * W[0]; 	     
		prob_sum += W[6];
	      };
	    };  
	    //   again time to decide what next

	    if(verboseLevel > 2){
	      G4cout << " wn " << W[0] << " wp " << W[1] << " wd " << W[2] << endl
		     << " wh3 " << W[3] << " wt " << W[4] << " whe4 " << W[5] << endl
		     << " wfi " << W[6] << G4endl;
	    }

	    G4int icase = -1;

	    if(prob_sum < prob_cut_off) { // photon emission chain

	      G4double UCR0 = 2.5 + 150.0 / A;
	      G4double T00 = 1.0 / (sqrt(u1 / UCR0) - 1.25 / UCR0);
	      G4int itry_gam = 0;

	      while(EEXS > cut_off_energy && try_again) {
		itry_gam++;

		G4int itry = 0;
		G4double T04 = 4.0 * T00;
		G4double FMAX;

		if(T04 < EEXS) {
		  FMAX = pow(T04, 4) * exp((EEXS - T04) / T00);	        
		}
		else {
		  FMAX = pow(EEXS, 4);
		}; 

		G4double S;

		while(itry < itry_max) {
		  itry++;
		  S = EEXS * inuclRndm();

		  G4double X1 = pow(S, 4) * exp((EEXS - S) / T00);

		  if(X1 > FMAX * inuclRndm()) break;
		};
		if(itry < itry_max) {
		  //      new photon escape

		  G4InuclElementaryParticle particle(10);
		  G4double pmod = 0.001 * S;
		  G4std::vector<G4double> mom(4);
		  G4std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
		  G4double FI = randomPHI();
		  G4double P1 = pmod * COS_SIN.second;

		  mom[1] = P1 * cos(FI);
		  mom[2] = P1 * sin(FI);
		  mom[3] = pmod * COS_SIN.first;
		  mom[0] = pmod;

		  G4std::vector<G4double> mom_at_rest(4);

		  for(G4int i = 1; i < 4; i++) mom_at_rest[i] = -mom[i];
		  mom_at_rest[0] = sqrt(mom_at_rest[1] * mom_at_rest[1] +
					mom_at_rest[2] * mom_at_rest[2] + 
					mom_at_rest[3] * mom_at_rest[3] +
					nuc_mass * nuc_mass); 

		  G4std::vector<G4double> part_mom = 
		    toTheNucleiSystemRestFrame.backToTheLab(mom);

		  part_mom[0] = sqrt(part_mom[1] * part_mom[1] +
				     part_mom[2] * part_mom[2] + 
				     part_mom[3] * part_mom[3]); 

		  G4std::vector<G4double> ex_mom = 
		    toTheNucleiSystemRestFrame.backToTheLab(mom_at_rest);

		  ex_mom[0] = sqrt(ex_mom[1] * ex_mom[1] + 
				   ex_mom[2] * ex_mom[2]
				   + ex_mom[3] * ex_mom[3] + 
				   nuc_mass * nuc_mass);			
		  EEXS_new = 1000.0 * (PEX[0] + 0.001 * EEXS - 
				       part_mom[0] - ex_mom[0]);
		  if(EEXS_new > 0.0) { // everything ok
		    PEX = ex_mom;
		    EEXS = EEXS_new;
		    particle.setMomentum(part_mom);
		    output.addOutgoingParticle(particle);
		    for(G4int i = 0; i < 4; i++) ppout[i] += part_mom[i];
		  }
		  else {
		    if(itry_gam == itry_gam_max) try_again = false;
		  };
		}
		else {
		  try_again = false;
		}; 
	      };
	      try_again = false;
	    }
	    else {

	      G4double SL = prob_sum * inuclRndm();
	      G4double S1 = 0.0;

	      for(G4int i = 0; i < 7; i++) {
		S1 += W[i]; 	
		if(SL <= S1) {
		  icase = i;

		  break;
		};
	      };
	      if(icase < 6) { // particle or light nuclei escape

		G4double uc = 2.0 * sqrt(u[icase] * TM[icase]);
		G4double ur = (uc > huge ? exp(huge) : exp(uc));
		G4double d1 = 1.0 / ur;
		G4double d2 = 1.0 / (ur - 1.0);	    
		G4int itry1 = 0;
		G4bool bad = true;

		while(itry1 < itry_max && bad) {
		  itry1++; 
	      
		  G4int itry = 0;
		  G4double EPR = -1.0;
		  G4double S = 0.0;

		  while(itry < itry_max && EPR < 0.0) {
		    itry++;

		    G4double uu = uc + log((1.0 - d1) * inuclRndm() + d2);

		    S = 0.5 * (uc * uc - uu * uu) / u[icase];
		    EPR = TM[icase] - S * A / (A - 1.0) + V[icase];
		  }; 
	    
		  if(EPR > 0.0 && S > 0.0) { // real escape
		    S = 0.001 * S;
		    if(icase < 2) { // particle escape

		      G4int ptype = 2 - icase;
		      G4InuclElementaryParticle particle(ptype);
		      G4double mass = particle.getMass();
		      //                       generate particle momentum
		      G4double pmod = sqrt((2.0 * mass + S) * S);
		      G4std::vector<G4double> mom(4);
		      G4std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
		      G4double FI = randomPHI();
		      G4double P1 = pmod * COS_SIN.second;

		      mom[1] = P1 * cos(FI);
		      mom[2] = P1 * sin(FI);
		      mom[3] = pmod * COS_SIN.first;

		      G4std::vector<G4double> mom_at_rest(4);

		      for(G4int i = 1; i < 4; i++) mom_at_rest[i] = -mom[i];

		      G4double new_nuc_mass = dummy_nuc.getNucleiMass(A1[icase],
								      Z1[icase]);

		      mom_at_rest[0] = sqrt(mom_at_rest[1] * mom_at_rest[1] +
					    mom_at_rest[2] * mom_at_rest[2] + 
					    mom_at_rest[3] * mom_at_rest[3] +
					    new_nuc_mass * new_nuc_mass); 
		      mom[0] = sqrt(mom[1] * mom[1] + mom[2] * mom[2] +
				    mom[3] * mom[3] + mass * mass);

		      G4std::vector<G4double> part_mom = 
		        toTheNucleiSystemRestFrame.backToTheLab(mom);

		      part_mom[0] = sqrt(part_mom[1] * part_mom[1] +
					 part_mom[2] * part_mom[2] + 
					 part_mom[3] * part_mom[3] +
					 mass * mass);

		      G4std::vector<G4double> ex_mom = 
		        toTheNucleiSystemRestFrame.backToTheLab(mom_at_rest);

		      ex_mom[0] = sqrt(ex_mom[1] * ex_mom[1] + 
				       ex_mom[2] * ex_mom[2] + 
				       ex_mom[3] * ex_mom[3] + 
				       new_nuc_mass * new_nuc_mass);			
		      EEXS_new = 1000.0 * (PEX[0] + 0.001 * EEXS - 
					   part_mom[0] - ex_mom[0]);

		      if(EEXS_new > 0.0) { // everything ok
			PEX = ex_mom;
			EEXS = EEXS_new;
			A = A1[icase];
			Z = Z1[icase]; 	      
			particle.setMomentum(part_mom);
			output.addOutgoingParticle(particle);
			for(G4int i = 0; i < 4; i++) ppout[i] += part_mom[i];
			bad = false;
		      };
		    }
		    else {

		      G4InuclNuclei nuclei(AN[icase], Q[icase]);
		      G4double mass = nuclei.getMass();
		      //                       generate particle momentum
		      G4double pmod = sqrt((2.0 * mass + S) * S);
		      G4std::vector<G4double> mom(4);
		      G4std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
		      G4double FI = randomPHI();
		      G4double P1 = pmod * COS_SIN.second;

		      mom[1] = P1 * cos(FI);
		      mom[2] = P1 * sin(FI);
		      mom[3] = pmod * COS_SIN.first;

		      G4std::vector<G4double> mom_at_rest(4);

		      for(G4int i = 1; i < 4; i++) mom_at_rest[i] = -mom[i];

		      G4double new_nuc_mass = dummy_nuc.getNucleiMass(A1[icase],
								      Z1[icase]);

		      mom_at_rest[0] = sqrt(mom_at_rest[1] * mom_at_rest[1] +
					    mom_at_rest[2] * mom_at_rest[2] + 
					    mom_at_rest[3] * mom_at_rest[3] +
					    new_nuc_mass * new_nuc_mass); 
		      mom[0] = sqrt(mom[1] * mom[1] + 
				    mom[2] * mom[2] +
				    mom[3] * mom[3] + 
				    mass * mass);

		      G4std::vector<G4double> part_mom = 
		        toTheNucleiSystemRestFrame.backToTheLab(mom);

		      part_mom[0] = sqrt(part_mom[1] * part_mom[1] +
					 part_mom[2] * part_mom[2] + 
					 part_mom[3] * part_mom[3] +
					 mass * mass);

		      G4std::vector<G4double> ex_mom = 
		        toTheNucleiSystemRestFrame.backToTheLab(mom_at_rest);

		      ex_mom[0] = sqrt(ex_mom[1] * ex_mom[1] + 
				       ex_mom[2] * ex_mom[2] +
				       ex_mom[3] * ex_mom[3] + 
				       new_nuc_mass * new_nuc_mass);			
		      EEXS_new = 1000.0 * (PEX[0] + 0.001 * EEXS - 
					   part_mom[0] - ex_mom[0]);
		      if(EEXS_new > 0.0) { // everything ok
			PEX = ex_mom;
			EEXS = EEXS_new;
			A = A1[icase];
			Z = Z1[icase]; 	      
			for(G4int i = 0; i < 4; i++) ppout[i] += part_mom[i];
			nuclei.setMomentum(part_mom);
			nuclei.setExitationEnergy(0.0);
			nuclei.setEnergy();
			output.addTargetFragment(nuclei);
			bad = false;
		      };
		    };
		  };
		};
		if(itry1 == itry_max || bad)  try_again = false;
	      }
	      else { // fission

		G4InuclNuclei nuclei(A, Z);        

		nuclei.setExitationEnergy(EEXS);

		if(verboseLevel > 2){
		  G4cout << " fission: A " << A << " Z " << Z << " eexs " << EEXS <<
		    " Wn " << W[0] << " Wf " << W[6] << G4endl;
		}

		G4CollisionOutput foutput = theFissioner->collide(0, &nuclei);
		G4std::vector<G4InuclNuclei> nuclea = foutput.getNucleiFragments();

		if(nuclea.size() == 2) { // fission o'k
		  //                convert back to the lab
		  for(G4int i = 0; i < 2; i++) {

		    G4std::vector<G4double> mom = nuclea[i].getMomentum();

		    mom = toTheNucleiSystemRestFrame.backToTheLab(mom);
		    nuclea[i].setMomentum(mom);
		    nuclea[i].setEnergy();
		  };	      

		  G4CollisionOutput output1 = this->collide(0, &nuclea[0]);

		  output.addOutgoingParticles(output1.getOutgoingParticles());
		  output.addTargetFragments(output1.getNucleiFragments());
		  output1 = this->collide(0, &nuclea[1]);
		  output.addOutgoingParticles(output1.getOutgoingParticles());
		  output.addTargetFragments(output1.getNucleiFragments());

		  return output;
		}
		else { // fission forbidden now
		  fission_open = false;
		}; 
	      }; 
	    }; 	     
	  }
	  else {
	    try_again = false;
	  }; 
	}; 
      };
      //   this time it's final nuclei
      if(itry_global == itry_global_max) G4cout << " ! itry_global " <<
					   itry_global_max << G4endl;

      G4std::vector<G4double> pnuc(4);

      for(G4int i = 1; i < 4; i++) pnuc[i] = pin[i] - ppout[i];

      G4InuclNuclei nuclei(pnuc, A, Z);

      nuclei.setEnergy();
      pnuc = nuclei.getMomentum(); 

      G4double eout = pnuc[0] + ppout[0];  
      G4double eex_real = 1000.0 * (pin[0] - eout);        

      nuclei.setExitationEnergy(eex_real);
      output.addTargetFragment(nuclei);
    };
  } 
  else {

    G4cout << " EquilibriumEvaporator -> target is not nuclei " << G4endl;    

  }; 

  return output;
}		     

G4bool G4EquilibriumEvaporator::timeToBigBang(G4double a, 
					      G4double z, 
					      G4double e) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::timeToBigBang" << G4endl;
  }

  const G4double be_cut = 3.0;

  G4bool bigb = true;

  if(a >= 12 && z >= 6.0 && z < 3.0 * (a - z)) {
    bigb = false;  
  }
  else {
    if(e < be_cut * bindingEnergy(a, z)) bigb = false;
  };

  return bigb;
}

G4bool G4EquilibriumEvaporator::goodRemnant(G4double a, 
					    G4double z) const {
  
  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::goodRemnant" << G4endl;
  }

  return a > 1.0 && z > 0.0 && a > z;
}

G4double G4EquilibriumEvaporator::getQF(G4double x, 
					G4double x2, 
					G4double a,
					G4double z, 
					G4double e) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::getQF" << G4endl;
  }
  
  const G4double QFREP[72] = {  
    //     TL201 *     *   *    *
    //      1    2     3   4    5
    22.5, 22.0, 21.0, 21.0, 20.0,
    //     BI209 BI207 PO210 AT213 *    TH234
    //      6     7    8     9     10   11
    20.6, 20.6, 18.6, 15.8, 13.5, 6.5,
    //     TH233 TH232 TH231 TH230 TX229 PA233 PA232 PA231 PA230 U240
    //     12    13    14    15    16    17    18    19    20    21
    6.65, 6.22, 6.27, 6.5,  6.7,  6.2,  6.25, 5.9,  6.1,  5.75,
    //     U239 U238 U237  U236 U235 U234 U233 U232 U231
    //     22   23   24    25   26   27   28   29   30
    6.46, 5.7, 6.28, 5.8, 6.15, 5.6, 5.8, 5.2, 5.8,
    //     NP238 NP237 NP236 NP235 PU245 NP234  PU244 NP233
    //     31    32    33    34    35    36     37    38
    6.2 , 5.9 , 5.9,  6.0,  5.8,  5.7,   5.4,  5.4,
    //     PU242 PU241 PU240 PU239 PU238 AM247 PU236 AM245 AM244 AM243
    //     39    40    41    42    43    44    45    46    47    48
    5.6,  6.1,  5.57, 6.3,  5.5,  5.8,  4.7,  6.2,  6.4,  6.2,
    //     AM242 AM241 AM240 CM250 AM239 CM249 CM248 CM247 CM246
    //     49    50    51    52    53    54    55    56    57
    6.5,  6.2,  6.5,  5.3,  6.4,  5.7,  5.7,  6.2,  5.7,
    //     CM245 CM244 CM243 CM242 CM241 BK250 CM240
    //     58    59    60    61    62    63    64
    6.3,  5.8,  6.7,  5.8,  6.6,  6.1,  4.3,
    //     BK249 CF252 CF250 CF248 CF246 ES254 ES253 FM254
    //     65    66    67    68    69    70    71    72
    6.2,  3.8,  5.6,  4.0,  4.0,  4.2,  4.2,  3.5 };
     
  const G4double XREP[72] = {
    //      1      2     3      4      5
    0.6761, 0.677, 0.6788, 0.6803, 0.685,
    //      6     7     8     9     10     11
    0.6889, 0.6914, 0.6991, 0.7068, 0.725, 0.7391,
    //     12  13    14   15   16    17  18    19    20    21
    0.74, 0.741, 0.742, 0.743, 0.744, 0.7509, 0.752, 0.7531, 0.7543, 0.7548,
    //     22    23    24
    0.7557, 0.7566, 0.7576,
    //      25     26   27    28    29   30   31    32     33    34
    0.7587, 0.7597, 0.7608, 0.762, 0.7632, 0.7644, 0.7675, 0.7686, 0.7697, 0.7709,
    //      35    36    37    38    39   40    41
    0.7714, 0.7721, 0.7723, 0.7733, 0.7743, 0.7753, 0.7764,
    //      42    43    44    45    46    47    48   49
    0.7775, 0.7786, 0.7801, 0.781, 0.7821, 0.7831, 0.7842, 0.7852,
    //     50     51    52    53    54    55    56    57    58
    0.7864, 0.7875, 0.7880, 0.7887, 0.7889, 0.7899, 0.7909, 0.7919, 0.7930,
    //      59    60    61    62    63    64
    0.7941, 0.7953, 0.7965, 0.7977, 0.7987, 0.7989,
    //      65    66    67    68    69    70    71    72
    0.7997, 0.8075, 0.8097, 0.8119, 0.8143, 0.8164, 0.8174, 0.8274 };

  const G4double G0 = 20.4;
  const G4double XMIN = 0.6761;
  const G4double XMAX = 0.8274;

  G4double QFF = 0.0;

  if(x < XMIN || x > XMAX) {

    G4double X1 = 1.0 - 0.02 * x2;
    G4double FX=(0.73 + (3.33 * X1 - 0.66) * X1) * pow(X1, 3);

    QFF = G0 * FX * pow(a, 2.0 / 3.0);
  }
  else {
    for(G4int i = 1; i < 72; i++) {
      if(x <= XREP[i]) {
	QFF = QFREP[i - 1] + (QFREP[i] - QFREP[i - 1]) * (x - XREP[i - 1])/
	  (XREP[i] - XREP[i - 1]);
	break;
      };
    };
  };
  if(QFF < 0.0) QFF = 0.0;

  return QFF; 
}

G4double G4EquilibriumEvaporator::getAF(G4double x, 
					G4double a, 
					G4double z, 
					G4double e) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::getAF" << G4endl;
  }

  // ugly parameterisation to fit the experimental fission cs for Hg - Bi nuclei
  G4double AF = 1.285 * (1.0 - e / 1100.0);

  if(AF < 1.06) AF = 1.06;
  // if(AF < 1.08) AF = 1.08;

  return AF;
}	

G4double G4EquilibriumEvaporator::getPARLEVDEN(G4double A, 
					       G4double Z) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::getPARLEVDEN" << G4endl;
  }

  const G4double par = 0.125;

  return par;
}

G4double G4EquilibriumEvaporator::getE0(G4double A) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::getE0" << G4endl;
  }

  const G4double e0 = 200.0;   

  return e0;   
}
