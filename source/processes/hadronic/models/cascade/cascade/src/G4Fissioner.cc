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
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100318  M. Kelsey -- Bug fix setting mass with G4LV
// 20100319  M. Kelsey -- Use new generateWithRandomAngles for theta,phi stuff;
//		eliminate some unnecessary std::pow()
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class
// 20100622  M. Kelsey -- Use local "bindingEnergy()" to call through
// 20100711  M. Kelsey -- Add energy-conservation checking, reduce if-cascades
// 20100713  M. Kelsey -- Don't add excitation energy to mass (already there)
// 20100714  M. Kelsey -- Move conservation checking to base class
// 20100728  M. Kelsey -- Make fixed arrays static, move G4FissionStore to data
//		member and reuse
// 20100914  M. Kelsey -- Migrate to integer A and Z
// 20110214  M. Kelsey -- Follow G4InuclParticle::Model enumerator migration
// 20110801  M. Kelsey -- Replace constant-size std::vector's w/C arrays
// 20110922  M. Kelsey -- Follow G4InuclParticle::print(ostream&) migration
// 20120517  A. Ribon  -- Removed static in local vectors for reproducibility
// 20130622  Inherit from G4CascadeDeexciteBase, move to deExcite() interface
//		with G4Fragment
// 20130628  Replace local list of fragments with use of output G4Fragments
// 20150608  M. Kelsey -- Label all while loops as terminating.
// 20150619  M. Kelsey -- Replace std::exp with G4Exp
// 20150622  M. Kelsey -- For new G4cbrt(int), move powers of A outside.

#include "G4Fissioner.hh"
#include "G4CollisionOutput.hh"
#include "G4Exp.hh"
#include "G4Fragment.hh"
#include "G4HadTmpUtil.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4FissionStore.hh"
#include "G4FissionConfiguration.hh"
#include "G4InuclSpecialFunctions.hh"

using namespace G4InuclSpecialFunctions;


void G4Fissioner::deExcite(const G4Fragment& target,
			   G4CollisionOutput& output) {
  if (verboseLevel) {
    G4cout << " >>> G4Fissioner::deExcite" << G4endl;
  }

  if (verboseLevel > 1) 
    G4cout << " Fissioner input\n" << target << G4endl;

  // Initialize buffer for fission possibilities
  fissionStore.setVerboseLevel(verboseLevel);
  fissionStore.clear();

  getTargetData(target);
  G4double A13 = G4cbrt(A);
  G4double mass_in = PEX.m();
  G4double e_in = mass_in; 		// Mass includes excitation
  G4double PARA = 0.055 * A13*A13 * (G4cbrt(A-Z) + G4cbrt(Z));
  G4double TEM = std::sqrt(EEXS / PARA);
  G4double TETA = 0.494 * A13 * TEM;
  
  TETA = TETA / std::sinh(TETA);
  
  if (A < 246) PARA += (nucleiLevelDensity(A) - PARA) * TETA;
  
  G4int A1 = A/2 + 1;
  G4int Z1;
  G4int A2 = A - A1;

  G4double ALMA = -1000.0;
  G4double DM1 = bindingEnergy(A,Z);
  G4double EVV = EEXS - DM1;
  G4double DM2 = bindingEnergyAsymptotic(A, Z);
  G4double DTEM = (A < 220 ? 0.5 : 1.15);
  
  TEM += DTEM;
  
  G4double AL1[2] = { -0.15, -0.15 };
  G4double BET1[2] = { 0.05, 0.05 };

  G4double R12 = G4cbrt(A1) + G4cbrt(A2); 
  
  for (G4int i = 0; i < 50 && A1 > 30; i++) {
    A1--;
    A2 = A - A1;
    G4double X3 = 1.0 / G4cbrt(A1);
    G4double X4 = 1.0 / G4cbrt(A2);
    Z1 = G4lrint(getZopt(A1, A2, Z, X3, X4, R12) - 1.);
    G4double EDEF1[2];
    G4int Z2 = Z - Z1;
    G4double VPOT, VCOUL;
    
    potentialMinimization(VPOT, EDEF1, VCOUL, A1, A2, Z1, Z2, AL1, BET1, R12);
    
    G4double DM3 = bindingEnergy(A1,Z1);
    G4double DM4 = bindingEnergyAsymptotic(A1, Z1);
    G4double DM5 = bindingEnergy(A2,Z2);
    G4double DM6 = bindingEnergyAsymptotic(A2, Z2);
    G4double DMT1 = DM4 + DM6 - DM2;
    G4double DMT = DM3 + DM5 - DM1;
    G4double EZL = EEXS + DMT - VPOT;
    
    if(EZL > 0.0) { // generate fluctuations
      //  faster, using randomGauss
      G4double C1 = std::sqrt(getC2(A1, A2, X3, X4, R12) / TEM);
      G4double DZ = randomGauss(C1);
      
      DZ = DZ > 0.0 ? DZ + 0.5 : -std::fabs(DZ - 0.5);
      Z1 += G4int(DZ);
      Z2 -= G4int(DZ);
      
      G4double DEfin = randomGauss(TEM);	
      G4double EZ = (DMT1 + (DMT - DMT1) * TETA - VPOT + DEfin) / TEM;
      
      if (EZ >= ALMA) ALMA = EZ;
      G4double EK = VCOUL + DEfin + 0.5 * TEM;
      G4double EV = EVV + bindingEnergy(A1,Z1) + bindingEnergy(A2,Z2) - EK;
      
      if (EV > 0.0) fissionStore.addConfig(A1, Z1, EZ, EK, EV);
    };
  };
  
  std::size_t store_size = fissionStore.size();
  if (store_size == 0) return;		// No fission products

  G4FissionConfiguration config = 
    fissionStore.generateConfiguration(ALMA, inuclRndm());
  
  A1 = G4int(config.afirst);
  A2 = A - A1;
  Z1 = G4int(config.zfirst);
  
  G4int Z2 = Z - Z1;
  
  G4double mass1 = G4InuclNuclei::getNucleiMass(A1,Z1);
  G4double mass2 = G4InuclNuclei::getNucleiMass(A2,Z2);
  G4double EK = config.ekin;
  G4double pmod = std::sqrt(0.001 * EK * mass1 * mass2 / mass_in);
  
  G4LorentzVector mom1 = generateWithRandomAngles(pmod, mass1);
  G4LorentzVector mom2; mom2.setVectM(-mom1.vect(), mass2);
  
  G4double e_out = mom1.e() + mom2.e();
  G4double EV = 1000.0 * (e_in - e_out) / A;
  if (EV <= 0.0) return;		// No fission energy

  G4double EEXS1 = EV*A1;
  G4double EEXS2 = EV*A2;

  // Pass only last two nuclear fragments
  output.addRecoilFragment(makeFragment(mom1, A1, Z1, EEXS1));
  output.addRecoilFragment(makeFragment(mom2, A2, Z2, EEXS2));
}

G4double G4Fissioner::getC2(G4int A1, 
			    G4int A2, 
			    G4double X3, 
			    G4double X4, 
			    G4double R12) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Fissioner::getC2" << G4endl;
  }

  G4double C2 = 124.57 * (1.0 / A1 + 1.0 / A2) + 0.78 * (X3 + X4) - 176.9 *
    ((X3*X3*X3*X3) + (X4*X4*X4*X4)) + 219.36 * (1.0 / (A1 * A1) + 1.0 / (A2 * A2)) - 1.108 / R12;

  return C2;   
}

G4double G4Fissioner::getZopt(G4int A1, 
			      G4int A2, 
			      G4int ZT, 
			      G4double X3, 
			      G4double X4, 
			      G4double R12) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Fissioner::getZopt" << G4endl;
  }

  G4double Zopt = (87.7 * (X4 - X3) * (1.0 - 1.25 * (X4 + X3)) +
		   ZT * ((124.57 / A2 + 0.78 * X4 - 176.9 * (X4*X4*X4*X4) + 219.36 / (A2 * A2)) - 0.554 / R12)) /
    getC2(A1, A2, X3, X4, R12);

  return Zopt;
}	     

void G4Fissioner::potentialMinimization(G4double& VP, 
					G4double( &ED)[2],
					G4double& VC, 
					G4int AF, 
					G4int AS, 
					G4int ZF, 
					G4int ZS,
					G4double AL1[2], 
					G4double BET1[2], 
					G4double& R12) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4Fissioner::potentialMinimization" << G4endl;
  }

  const G4double huge_num = 2.0e35;
  const G4int itry_max = 2000;
  const G4double DSOL1 = 1.0e-6;
  const G4double DS1 = 0.3;
  const G4double DS2 = 1.0 / DS1 / DS1; 
  G4int A1[2] = { AF, AS };
  G4int Z1[2] = { ZF, ZS };
  G4double D = 1.01844 * ZF * ZS;
  G4double D0 = 1.0e-3 * D;
  G4double R[2];
  R12 = 0.0;
  G4double C[2];
  G4double F[2];
  G4double Y1;
  G4double Y2;
  G4int i;

  for (i = 0; i < 2; i++) {
    R[i] = G4cbrt(A1[i]);
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
  G4double AA[4][4];
  G4double B[4];
  G4int itry = 0;

  while (itry < itry_max) {	/* Loop checking 08.06.2015 MHK */
    itry++;
    G4double S = 0.0;

    for (i = 0; i < 2; i++) {
      S += R[i] * (1.0 + AL1[i] + BET1[i] - 0.257 * AL1[i] * BET1[i]);
    };
    R12 = 0.0;
    Y1 = 0.0;
    Y2 = 0.0;

    for (i = 0; i < 2; i++) {
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
 
    for (i = 0; i < 2; i++) {
      RAL[i] = -R[i] * (1.0 - 0.6 * X[i]) + SAL[i] * Y3;
      RBE[i] =  R[i] * (1.0 - 0.429 * X1[i]) + SBE[i] * Y3;
    };

    G4double DX1;
    G4double DX2;
  
    for (i = 0; i < 2; i++) {

      for (G4int j = 0; j < 2; j++) {
	G4double DEL1 = i == j ? 1.0 : 0.0;
	DX1 = 0.0;
	DX2 = 0.0;

	if (std::fabs(AL1[i]) >= DS1) {
	  G4double XXX = AL1[i] * AL1[i] * DS2;
	  G4double DEX = XXX > 100.0 ? huge_num : G4Exp(XXX);
	  DX1 = 2.0 * (1.0 + 2.0 * AL1[i] * AL1[i] * DS2) * DEX * DS2;
	};

	if (std::fabs(BET1[i]) >= DS1) {
	  G4double XXX = BET1[i] * BET1[i] * DS2;
	  G4double DEX = XXX > 100.0 ? huge_num : G4Exp(XXX);
	  DX2 = 2.0 * (1.+2.0 * BET1[i] * BET1[i] * DS2) * DEX * DS2;
	};

	G4double DEL = 2.0e-3 * DEL1;
	AA[i][j] = R3 * RBE[i] * RBE[j] - 
	  R2 * (-0.6 * 
		(X1[i] * SAL[j] + 
		 X1[j] * SAL[i]) + SAL[i] * SAL[j] * Y4) + 
	  DEL * C[i] + DEL1 * DX1;
	G4int i1 = i + 2;
	G4int j1 = j + 2;
	AA[i1][j1] = R3 * RBE[i] * RBE[j] 
	  - R2 * (0.857 * 
		  (X2[i] * SBE[j] + 
		   X2[j] * SBE[i]) + SBE[i] * SBE[j] * Y4) +
	  DEL * F[i] + DEL1 * DX2;
	AA[i][j1] = R3 * RAL[i] * RBE[j] - 
	  R2 * (0.857 * 
		(X2[j] * SAL[i] - 
		 0.6 * X1[i] * SBE[j]) + SBE[j] * SAL[i] * Y4 - 
		0.257 * R[i] * Y3 * DEL1);
	AA[j1][i] = AA[i][j1];   	 
      };
    };

    for (i = 0; i < 2; i++) {
      DX1 = 0.0;
      DX2 = 0.0;

      if (std::fabs(AL1[i]) >= DS1) DX1 = 2.0 * AL1[i] * DS2 * G4Exp(AL1[i] * AL1[i] * DS2);

      if (std::fabs(BET1[i]) >= DS1) DX2 = 2.0 * BET1[i] * DS2 * G4Exp(BET1[i] * BET1[i] * DS2);
      B[i] =     R2 * RAL[i] - 2.0e-3 * C[i] * AL1[i] + DX1;
      B[i + 2] = R2 * RBE[i] - 2.0e-3 * F[i] * BET1[i] + DX2;
    };

    G4double ST = 0.0;
    G4double ST1 = 0.0;

    for (i = 0; i < 4; i++) {
      ST += B[i] * B[i];

      for (G4int j = 0; j < 4; j++) ST1 += AA[i][j] * B[i] * B[j];
    };

    G4double STEP = ST / ST1;
    G4double DSOL = 0.0;

    for (i = 0; i < 2; i++) {
      AL1[i] += B[i] * STEP;
      BET1[i] += B[i + 2] * STEP;
      DSOL += B[i] * B[i] + B[i + 2] * B[i + 2]; 
    };
    DSOL = std::sqrt(DSOL);

    if (DSOL < DSOL1) break;
  };

  if (verboseLevel > 3) {
  if (itry == itry_max) 
    G4cout << " maximal number of iterations in potentialMinimization " << G4endl
	   << " A1 " << AF << " Z1 " << ZF << G4endl; 

  };

  for (i = 0; i < 2; i++) ED[i] = F[i] * BET1[i] * BET1[i] + C[i] * AL1[i] * AL1[i]; 

  VC = D / R12;
  VP = VC + ED[0] + ED[1];
}
