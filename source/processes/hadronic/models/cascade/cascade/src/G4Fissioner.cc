#include "G4Fissioner.hh"
#include "G4InuclNuclei.hh"
#include "G4FissionStore.hh"
#include "G4FissionConfiguration.hh"

G4Fissioner::G4Fissioner()
  : verboseLevel(1) {
  
  if (verboseLevel > 3) {
    G4cout << " >>> G4Fissioner::G4Fissioner" << G4endl;
  }
}

G4CollisionOutput G4Fissioner::collide(G4InuclParticle* bullet,
				       G4InuclParticle* target) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Fissioner::collide" << G4endl;
  }

  const G4double one_third = 1.0/3.0;
  const G4double two_thirds = 2.0/3.0;
  //  const G4int itry_max = 1000;

  G4CollisionOutput output;

  if(G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target)) {

    if (verboseLevel > 3) {
      G4cout << " Fissioner input " << G4endl;

      nuclei_target->printParticle();
    }

    G4double A = nuclei_target->getA();
    G4double Z = nuclei_target->getZ();
    G4double EEXS = nuclei_target->getExitationEnergy();
    G4double mass_in = nuclei_target->getMass();
    G4double e_in = mass_in + 0.001 * EEXS;
    G4double PARA = 0.055 * pow(A, two_thirds) * (pow(A - Z, one_third) + pow(Z, one_third));
    G4double TEM = sqrt(EEXS / PARA);
    G4double TETA = 0.494 * pow(A, one_third) * TEM;

    TETA = TETA / sinh(TETA);

    if(A < 246.0) PARA += (nucleiLevelDensity(A) - PARA) * TETA;

    G4double A1 = int(A / 2.0 + 1.1);
    G4double Z1;
    G4double A2 = A - A1;
    G4double ALMA = -1000.0;
    G4double DM1 = bindingEnergy(A, Z);
    G4double EVV = EEXS - DM1;
    G4double DM2 = bindingEnergyAsymptotic(A, Z);
    G4double DTEM = (A < 220.0 ? 0.5 : 1.15);

    TEM += DTEM;
  
    G4std::vector<G4double> AL1(2, -0.15);
    G4std::vector<G4double> BET1(2, 0.05);
    G4FissionStore fissionStore;
    G4double R12 = pow(A1, one_third) + pow(A2, one_third); 
  
    for(G4int i = 0; i < 50 && A1 > 30.0; i++) {
      A1 -= 1.0;
      A2 = A - A1;

      G4double X3 = 1.0 / pow(A1, one_third);
      G4double X4 = 1.0 / pow(A2, one_third);

      Z1 = int(getZopt(A1, A2, Z, X3, X4, R12)) - 1.0;

      G4std::vector<G4double> EDEF1(2);
      G4double Z2 = Z - Z1;
      G4double VPOT, VCOUL;

      potentialMinimization(VPOT, EDEF1, VCOUL, A1, A2, Z1, Z2, AL1, BET1, R12);

      G4double DM3 = bindingEnergy(A1, Z1);
      G4double DM4 = bindingEnergyAsymptotic(A1, Z1);
      G4double DM5 = bindingEnergy(A2, Z2);
      G4double DM6 = bindingEnergyAsymptotic(A2, Z2);
      G4double DMT1 = DM4 + DM6 - DM2;
      G4double DMT = DM3 + DM5 - DM1;
      G4double EZL = EEXS + DMT - VPOT;
    
      if(EZL > 0.0) { // generate fluctuations
	//  faster, using randomGauss

	G4double C1 = sqrt(getC2(A1, A2, X3, X4, R12) / TEM);
	G4double DZ = randomGauss(C1);

	DZ = DZ > 0.0 ? int(DZ + 0.5) : -int(fabs(DZ - 0.5));
	Z1 += DZ;
	Z2 -= DZ;

	G4double DEfin = randomGauss(TEM);	
	G4double EZ = (DMT1 + (DMT - DMT1) * TETA - VPOT + DEfin) / TEM;

	if(EZ >= ALMA) ALMA = EZ;

	G4double EK = VCOUL + DEfin + 0.5 * TEM;
	G4double EV = EVV + bindingEnergy(A1, Z1) + bindingEnergy(A2, Z2) - EK;
       
	if(EV > 0.0) fissionStore.addConfig(A1, Z1, EZ, EK, EV);
      };
    };

    G4int store_size = fissionStore.size();

    if(store_size > 0) {

      G4FissionConfiguration config = 
	fissionStore.generateConfiguration(ALMA, inuclRndm());

      A1 = config.afirst;
      A2 = A - A1;
      Z1 = config.zfirst;

      G4double Z2 = Z - Z1;
      G4InuclNuclei nuclei1(A1, Z1);
      G4InuclNuclei nuclei2(A2, Z2);        
      G4double mass1 = nuclei1.getMass();
      G4double mass2 = nuclei2.getMass();
      G4double EK = config.ekin;
      G4double pmod = sqrt(0.001 * EK * mass1 * mass2 / mass_in);
      pair<G4double, G4double> COS_SIN = randomCOS_SIN();
      G4double Fi = randomPHI();
      G4double P1 = pmod * COS_SIN.second;
      G4std::vector<G4double> mom1(4);
      G4std::vector<G4double> mom2(4);

      mom1[1] = P1 * cos(Fi);
      mom1[2] = P1 * sin(Fi);
      mom1[3] = pmod * COS_SIN.first;

      for(G4int i = 1; i < 4; i++) mom2[i] = -mom1[i];

      G4double e_out = sqrt(pmod * pmod + mass1 * mass1) + 
	sqrt(pmod * pmod + mass2 * mass2);
      G4double EV = 1000.0 * (e_in - e_out) / A;

      if(EV > 0.0) {

	G4double EEXS1 = EV*A1;
	G4double EEXS2 = EV*A2;
	G4InuclNuclei nuclei1(mom1, A1, Z1);        

	nuclei1.setExitationEnergy(EEXS1);
	nuclei1.setEnergy();
	output.addTargetFragment(nuclei1);

	G4InuclNuclei nuclei2(mom2, A2, Z2);        

	nuclei2.setExitationEnergy(EEXS2);
	nuclei2.setEnergy();
	output.addTargetFragment(nuclei2);

	if (verboseLevel > 3) {
	  nuclei1.printParticle();
	  nuclei2.printParticle();
	}
      };
    };
  }
  else {

    G4cout << " Fissioner -> target is not nuclei " << G4endl;    

  }; 

  return output;
}

G4double G4Fissioner::getC2(G4double A1, 
			    G4double A2, 
			    G4double X3, 
			    G4double X4, 
			    G4double R12) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Fissioner::getC2" << G4endl;
  }

  G4double C2 = 124.57 * (1.0 / A1 + 1.0 / A2) + 0.78 * (X3 + X4) - 176.9 *
    (pow(X3, 4) + pow(X4, 4)) + 219.36 * (1.0 / (A1 * A1) + 1.0 / (A2 * A2)) - 1.108 / R12;

  return C2;   
}

G4double G4Fissioner::getZopt(G4double A1, 
			      G4double A2, 
			      G4double ZT, 
			      G4double X3, 
			      G4double X4, 
			      G4double R12) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Fissioner::getZopt" << G4endl;
  }

  G4double Zopt = (87.7 * (X4 - X3) * (1.0 - 1.25 * (X4 + X3)) +
		   ZT * ((124.57 / A2 + 0.78 * X4 - 176.9 * pow(X4, 4) + 219.36 / (A2 * A2)) - 0.554 / R12)) /
    getC2(A1, A2, X3, X4, R12);

  return Zopt;
}	     

void G4Fissioner::potentialMinimization(G4double& VP, 
					G4std::vector<G4double> & ED,
					G4double& VC, 
					G4double AF, 
					G4double AS, 
					G4double ZF, 
					G4double ZS,
					G4std::vector<G4double>& AL1, 
					G4std::vector<G4double>& BET1, 
					G4double& R12) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Fissioner::potentialMinimization" << G4endl;
  }

  const G4double huge = 2.0e35;
  const G4double one_third = 1.0 / 3.0;
  //  const G4double two_thirds = 2.0 / 3.0;
  const G4int itry_max = 2000;
  const G4double DSOL1 = 1.0e-6;
  const G4double DS1 = 0.3;
  const G4double DS2 = 1.0 / DS1 / DS1; 

  G4double A1[2];

  A1[0] = AF;
  A1[1] = AS;

  G4double Z1[2];

  Z1[0] = ZF;
  Z1[1] = ZS;

  G4double D = 1.01844 * ZF * ZS;
  G4double D0 = 1.0e-3 * D;
  G4double R[2];

  R12 = 0.0;

  G4double C[2];
  G4double F[2];
  G4double Y1;
  G4double Y2;

  for(G4int i = 0; i < 2; i++) {
    R[i] = pow(A1[i], one_third);
    Y1 = R[i] * R[i];
    Y2 = Z1[i] * Z1[i] / R[i];
    C[i] = 6.8 * Y1 - 0.142 * Y2;
    F[i] = 12.138 * Y1 - 0.145 * Y2; 
  };

  G4double SAL[2];
  G4double SBE[2];
  G4double X[2];
  G4double X1[2];
  G4double X2[2];
  G4double RAL[2];
  G4double RBE[2];
  G4double A[4][4];
  G4double B[4];
  G4int itry = 0;

  while(itry < itry_max) {
    itry++;

    G4double S = 0.0;

    for(G4int i = 0; i < 2; i++) {
      S += R[i] * (1.0 + AL1[i] + BET1[i] - 0.257 * AL1[i] * BET1[i]);
    };
    R12 = 0.0;
    Y1 = 0.0;
    Y2 = 0.0;
    for(i = 0; i < 2; i++) {
      SAL[i] = R[i] * (1.0-0.257 * BET1[i]);
      SBE[i] = R[i] * (1.0-0.257 * AL1[i]);
      X[i] = R[i] / S;
      X1[i] = X[i] * X[i];
      X2[i] = X[i] * X1[i];
      Y1 += AL1[i] * X1[i];
      Y2 += BET1[i] * X2[i];
      R12 += R[i] * (1.0 - AL1[i] * (1.0 - 0.6 * X[i]) + BET1[i] * (1.0 - 0.429 * X1[i]));
    }; 

    G4double Y3 = -0.6 * Y1 + 0.857 * Y2;
    G4double Y4 = (1.2 * Y1 - 2.571 * Y2) / S;
    G4double R2 = D0 / (R12 * R12);
    G4double R3 = 2.0 * R2 / R12;
 
    for(i = 0; i < 2; i++) {
      RAL[i] = -R[i] * (1.0 - 0.6 * X[i]) + SAL[i] * Y3;
      RBE[i] =  R[i] * (1.0 - 0.429 * X1[i]) + SBE[i] * Y3;
    };

    G4double DX1;
    G4double DX2;
  
    for(i = 0; i < 2; i++) {
      for(G4int j = 0; j < 2; j++) {

	G4double DEL1 = i == j ? 1.0 : 0.0;

	DX1 = 0.0;
	DX2 = 0.0;
	if(fabs(AL1[i]) >= DS1) {

	  G4double XXX = AL1[i] * AL1[i] * DS2;
	  G4double DEX = XXX > 100.0 ? huge : exp(XXX);

	  DX1 = 2.0 * (1.0 + 2.0 * AL1[i] * AL1[i] * DS2) * DEX * DS2;
	};
	if(fabs(BET1[i]) >= DS1) {

	  G4double XXX = BET1[i] * BET1[i] * DS2;
	  G4double DEX = XXX > 100.0 ? huge : exp(XXX);

	  DX2 = 2.0 * (1.+2.0 * BET1[i] * BET1[i] * DS2) * DEX * DS2;
	};

	G4double DEL = 2.0e-3 * DEL1;

	A[i][j] = R3 * RBE[i] * RBE[j] - 
	  R2 * (-0.6 * 
		(X1[i] * SAL[j] + 
		 X1[j] * SAL[i]) + SAL[i] * SAL[j] * Y4) + 
	  DEL * C[i] + DEL1 * DX1;

	G4int i1 = i + 2;
	G4int j1 = j + 2;

	A[i1][j1] = R3 * RBE[i] * RBE[j] 
	  - R2 * (0.857 * 
		  (X2[i] * SBE[j] + 
		   X2[j] * SBE[i]) + SBE[i] * SBE[j] * Y4) +
	  DEL * F[i] + DEL1 * DX2;
	A[i][j1] = R3 * RAL[i] * RBE[j] - 
	  R2 * (0.857 * 
		(X2[j] * SAL[i] - 
		 0.6 * X1[i] * SBE[j]) + SBE[j] * SAL[i] * Y4 - 
		0.257 * R[i] * Y3 * DEL1);
	A[j1][i] = A[i][j1];   	 
      };
    };
    for(i = 0; i < 2; i++) {
      DX1 = 0.0;
      DX2 = 0.0;
      if(fabs(AL1[i]) >= DS1) DX1 = 2.0 * AL1[i] * DS2 * exp(AL1[i] * AL1[i] * DS2);
      if(fabs(BET1[i]) >= DS1) DX2 = 2.0 * BET1[i] * DS2 * exp(BET1[i] * BET1[i] * DS2);
      B[i] =     R2 * RAL[i] - 2.0e-3 * C[i] * AL1[i] + DX1;
      B[i + 2] = R2 * RBE[i] - 2.0e-3 * F[i] * BET1[i] + DX2;
    };

    G4double ST = 0.0;
    G4double ST1 = 0.0;

    for(i = 0; i < 4; i++) {
      ST += B[i] * B[i];
      for(G4int j = 0; j < 4; j++) ST1 += A[i][j] * B[i] * B[j];
    };

    G4double STEP = ST / ST1;
    G4double DSOL = 0.0;

    for(i = 0; i < 2; i++) {
      AL1[i] += B[i] * STEP;
      BET1[i] += B[i + 2] * STEP;
      DSOL += B[i] * B[i] + B[i + 2] * B[i + 2]; 
    };
    DSOL = sqrt(DSOL);
    if(DSOL < DSOL1) break;
  };
  if(itry == itry_max) 

    G4cout << " maximal number of iterations in potentialMinimization " << G4endl
	   << " A1 " << AF << " Z1 " << ZF << G4endl; 

  for(i = 0; i < 2; i++) 
    ED[i] = F[i] * BET1[i] * BET1[i] + C[i] * AL1[i] * AL1[i]; 
  VC = D / R12;
  VP = VC + ED[0] + ED[1];

}
