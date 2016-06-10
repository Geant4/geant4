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
// $Id: G4NonEquilibriumEvaporator.cc 71942 2013-06-28 19:08:11Z mkelsey $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100309  M. Kelsey -- Use new generateWithRandomAngles for theta,phi stuff;
//		eliminate some unnecessary std::pow()
// 20100412  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100413  M. Kelsey -- Pass buffers to paraMaker[Truncated]
// 20100517  M. Kelsey -- Inherit from common base class
// 20100617  M. Kelsey -- Remove "RUN" preprocessor flag and all "#else" code
// 20100622  M. Kelsey -- Use local "bindingEnergy()" function to call through.
// 20100701  M. Kelsey -- Don't need to add excitation to nuclear mass; compute
//		new excitation energies properly (mass differences)
// 20100713  M. Kelsey -- Add conservation checking, diagnostic messages.
// 20100714  M. Kelsey -- Move conservation checking to base class
// 20100719  M. Kelsey -- Simplify EEXS calculations with particle evaporation.
// 20100724  M. Kelsey -- Replace std::vector<> D with trivial D[3] array.
// 20100914  M. Kelsey -- Migrate to integer A and Z: this involves replacing
//		a number of G4double terms with G4int, with consequent casts.
// 20110214  M. Kelsey -- Follow G4InuclParticle::Model enumerator migration
// 20110922  M. Kelsey -- Follow G4InuclParticle::print(ostream&) migration
// 20120608  M. Kelsey -- Fix variable-name "shadowing" compiler warnings.
// 20121009  M. Kelsey -- Add some high-verbosity debugging output
// 20130622  Inherit from G4CascadeDeexciteBase, move to deExcite() interface
//		with G4Fragment
// 20130808  M. Kelsey -- Use new object-version of paraMaker, for thread safety
// 20130924  M. Kelsey -- Replace std::pow with G4Pow::powN() for CPU speed
// 20150608  M. Kelsey -- Label all while loops as terminating.

#include <cmath>

#include "G4NonEquilibriumEvaporator.hh"
#include "G4SystemOfUnits.hh"
#include "G4CollisionOutput.hh"
#include "G4Fragment.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4LorentzConvertor.hh"
#include "G4Pow.hh"

using namespace G4InuclSpecialFunctions;


G4NonEquilibriumEvaporator::G4NonEquilibriumEvaporator()
  : G4CascadeDeexciteBase("G4NonEquilibriumEvaporator"),
    theParaMaker(verboseLevel), theG4Pow(G4Pow::GetInstance()) {}


void G4NonEquilibriumEvaporator::deExcite(const G4Fragment& target,
					  G4CollisionOutput& output) {
  if (verboseLevel) {
    G4cout << " >>> G4NonEquilibriumEvaporator::deExcite" << G4endl;
  }

  if (verboseLevel>1) G4cout << " evaporating target:\n" << target << G4endl;
  
  const G4int a_cut = 5;
  const G4int z_cut = 3;

  const G4double eexs_cut = 0.1;

  const G4double coul_coeff = 1.4;
  const G4int itry_max = 1000;
  const G4double small_ekin = 1.0e-6;
  const G4double width_cut = 0.005;

  getTargetData(target);
  G4LorentzVector pin = PEX;		// Save original four-vector for later
  
  G4ExitonConfiguration config(target);  
  G4int QPP = config.protonQuasiParticles;
  G4int QNP = config.neutronQuasiParticles; 
  G4int QPH = config.protonHoles;
  G4int QNH = config.neutronHoles; 
  
  G4int QP = QPP + QNP;
  G4int QH = QPH + QNH;
  G4int QEX = QP + QH;
  
  G4InuclElementaryParticle dummy(small_ekin, 1);
  G4LorentzConvertor toTheExitonSystemRestFrame;
  //*** toTheExitonSystemRestFrame.setVerbose(verboseLevel);
  toTheExitonSystemRestFrame.setBullet(dummy);
  
  G4double EFN = FermiEnergy(A, Z, 0);
  G4double EFP = FermiEnergy(A, Z, 1);
  
  G4int AR = A - QP;
  G4int ZR = Z - QPP;  
  G4int NEX = QEX;
  G4LorentzVector ppout;
  G4bool try_again = (NEX > 0);
  
  // Buffer for parameter sets
  std::pair<G4double, G4double> parms;
  
  while (try_again) {			/* Loop checking 08.06.2015 MHK */
    if (A >= a_cut && Z >= z_cut && EEXS > eexs_cut) { // ok
      // update exiton system (include excitation energy!)
      G4double nuc_mass = G4InuclNuclei::getNucleiMass(A, Z, EEXS); 
      PEX.setVectM(PEX.vect(), nuc_mass);
      toTheExitonSystemRestFrame.setTarget(PEX);
      toTheExitonSystemRestFrame.toTheTargetRestFrame();
      
      if (verboseLevel > 2) {
	G4cout << " A " << A << " Z " << Z << " mass " << nuc_mass
	       << " EEXS " << EEXS << G4endl; 
      }
      
      G4double MEL = getMatrixElement(A);
      G4double E0 = getE0(A);
      G4double PL = getParLev(A, Z);
      G4double parlev = PL / A;
      G4double EG = PL * EEXS;
      
      if (QEX < std::sqrt(2.0 * EG)) { // ok
	if (verboseLevel > 3)
	  G4cout << " QEX " << QEX << " < sqrt(2*EG) " << std::sqrt(2.*EG)
		 << " NEX " << NEX << G4endl;
	
	theParaMaker.getTruncated(Z, parms);
	const G4double& AK1 = parms.first;
	const G4double& CPA1 = parms.second;
	
	G4double VP = coul_coeff * Z * AK1 / (G4cbrt(A-1) + 1.0) /
	  (1.0 + EEXS / E0);
	G4double DM1 = bindingEnergy(A,Z);
	G4double BN = DM1 - bindingEnergy(A-1,Z);
	G4double BP = DM1 - bindingEnergy(A-1,Z-1);
	G4double EMN = EEXS - BN;
	G4double EMP = EEXS - BP - VP * A / (A-1);
	G4double ESP = 0.0;

	if (verboseLevel > 3) {
	  G4cout << " AK1 " << AK1 << " CPA1 " << " VP " << VP
		 << "\n bind(A,Z) " << DM1 << " dBind(N) " << BN 
		 << " dBind(P) " << BP
		 << "\n EMN " << EMN << " EMP " << EMP << G4endl;
	}

	if (EMN > eexs_cut) { // ok
	  G4int icase = 0;
	  
	  if (NEX > 1) {
	    G4double APH = 0.25 * (QP * QP + QH * QH + QP - 3 * QH);
	    G4double APH1 = APH + 0.5 * (QP + QH);
	    ESP = EEXS / QEX;
	    G4double MELE = MEL / ESP / (A*A*A);

	    if (verboseLevel > 3)
	      G4cout << " APH " << APH << " APH1 " << APH1 << " ESP " << ESP
		     << G4endl;

	    if (ESP > 15.0) {
	      MELE *= std::sqrt(15.0 / ESP);
	    } else if(ESP < 7.0) {
	      MELE *= std::sqrt(ESP / 7.0);
	      if (ESP < 2.0) MELE *= std::sqrt(ESP / 2.0);
	    };    

	    G4double F1 = EG - APH;
	    G4double F2 = EG - APH1;

	    if (verboseLevel > 3)
	      G4cout << " MELE " << MELE << " F1 " << F1 << " F2 " << F2
		     << G4endl;
	    
	    if (F1 > 0.0 && F2 > 0.0) {
	      G4double F = F2 / F1;
	      G4double M1 = 2.77 * MELE * PL;
	      G4double D[3] = { 0., 0., 0. };
	      D[0] = M1 * F2 * F2 * theG4Pow->powN(F, NEX-1) / (QEX+1);
	      if (verboseLevel > 3) {
		G4cout << " D[0] " << D[0] << " with F " << F
		       << " powN(F,NEX-1) " << theG4Pow->powN(F, NEX-1)
		       << G4endl;
	      }

	      if (D[0] > 0.0) {
		
		if (NEX >= 2) {
		  D[1] = 0.0462 / parlev / G4cbrt(A) * QP * EEXS / QEX;
		  
		  if (EMP > eexs_cut) 
		    D[2] = D[1] * theG4Pow->powN(EMP/EEXS, NEX) * (1.0 + CPA1);
		  D[1] *= theG4Pow->powN(EMN/EEXS, NEX) * getAL(A);   

		  if (verboseLevel > 3) {
		    G4cout << " D[1] " << D[1] << " with powN(EMN/EEXS, NEX) "
			   << theG4Pow->powN(EMN/EEXS, NEX) << G4endl
			   << " D[2] " << D[2] << " with powN(EMP/EEXS, NEX) "
			   << theG4Pow->powN(EMP/EEXS, NEX) << G4endl;
		  }

		  if (QNP < 1) D[1] = 0.0;
		  if (QPP < 1) D[2] = 0.0;
		  
		  try_again = NEX > 1 && (D[1] > width_cut * D[0] || 
					  D[2] > width_cut * D[0]);
		  
		  if (try_again) {
		    G4double D5 = D[0] + D[1] + D[2];
		    G4double SL = D5 * inuclRndm();
		    G4double S1 = 0.;

		    if (verboseLevel > 3)
		      G4cout << " D5 " << D5 << " SL " << SL << G4endl;

		    for (G4int i = 0; i < 3; i++) {
		      S1 += D[i]; 	
		      if (SL <= S1) {
			icase = i;
			break;
		      }
		    }

		    if (verboseLevel > 3)
		      G4cout << " got icase " << icase << G4endl;
		  }				// if (try_again)
		}				// if (NEX >= 2)
	      } else try_again = false;		// if (D[0] > 0)
	    } else try_again = false;		// if (F1>0 && F2>0)
	  }					// if (NEX > 1)
	  
	  if (try_again) {
	    if (icase > 0) { 			// N -> N-1 with particle escape
	      if (verboseLevel > 3)
		G4cout << " try_again icase " << icase << G4endl;
	      
	      G4double V = 0.0;
	      G4int ptype = 0;
	      G4double B = 0.0;
	      
	      if (A < 3.0) try_again = false;
	      
	      if (try_again) { 
		
		if (icase == 1) { // neutron escape
		  if (verboseLevel > 3)
		    G4cout << " trying neutron escape" << G4endl;

		  if (QNP < 1) icase = 0;
		  else {
		    B = BN;
		    V = 0.0;
		    ptype = 2;		  
		  };    
		} else { // proton esape
		  if (verboseLevel > 3)
		    G4cout << " trying proton escape" << G4endl;

		  if (QPP < 1) icase = 0;
		  else {
		    B = BP;
		    V = VP;
		    ptype = 1;
		    
		    if (Z-1 < 1) try_again = false;
		  };   
		};
	        
		if (try_again && icase != 0) {
		  if (verboseLevel > 3)
		    G4cout << " ptype " << ptype << " B " << B << " V " << V
			   << G4endl;

		  G4double EB = EEXS - B;
		  G4double E = EB - V * A / (A-1);
		  
		  if (E < 0.0) icase = 0;
		  else {
		    G4double E1 = EB - V;
		    G4double EEXS_new = -1.;
		    G4double EPART = 0.0;
		    G4int itry1 = 0;
		    G4bool bad = true;
		    
		    /* Loop checking 08.06.2015 MHK */
		    while (itry1 < itry_max && icase > 0 && bad) {
		      itry1++;
		      G4int itry = 0;		    
		      
		      /* Loop checking 08.06.2015 MHK */
		      while (EEXS_new < 0.0 && itry < itry_max) {
			itry++;
			G4double R = inuclRndm();
			G4double X;
			
			if (NEX == 2) {
			  X = 1.0 - std::sqrt(R);
			  
			} else {		         
			  G4double QEX2 = 1.0 / QEX;
			  G4double QEX1 = 1.0 / (QEX-1);
			  X = theG4Pow->powA(0.5*R, QEX2);
			  if (verboseLevel > 3) {
			    G4cout << " R " << R << " QEX2 " << QEX2
				   << " powA(R, QEX2) " << X << G4endl;
			  }
			  
			  for (G4int i = 0; i < 1000; i++) {
			    G4double DX = X * QEX1 * 
			      (1.0 + QEX2 * X * (1.0 - R / theG4Pow->powN(X, NEX)) / (1.0 - X));
			    if (verboseLevel > 3) {
			      G4cout << " NEX " << NEX << " powN(X, NEX) "
				     << theG4Pow->powN(X, NEX) << G4endl;
			    }
			    
			    X -= DX;
			    
			    if (std::fabs(DX / X) < 0.01) break;  
			    
			  };
			}; 
			EPART = EB - X * E1;
			EEXS_new = EB - EPART * A / (A-1);
		      }	// while (EEXS_new < 0.0...
		      
		      if (itry == itry_max || EEXS_new < 0.0) {
			icase = 0;
			continue;
		      }
		      
		      if (verboseLevel > 2)
			G4cout << " particle " << ptype << " escape " << G4endl;
		      
		      EPART /= GeV; 		// From MeV to GeV
		      
		      G4InuclElementaryParticle particle(ptype);
		      particle.setModel(G4InuclParticle::NonEquilib);
		      
		      // generate particle momentum
		      G4double mass = particle.getMass();
		      G4double pmod = std::sqrt(EPART * (2.0 * mass + EPART));
		      G4LorentzVector mom = generateWithRandomAngles(pmod,mass);
		      
		      // Push evaporated paricle into current rest frame
		      mom = toTheExitonSystemRestFrame.backToTheLab(mom);
		      
		      // Adjust quasiparticle and nucleon counts
		      G4int QPP_new = QPP;
		      G4int QNP_new = QNP;
		      
		      G4int A_new = A-1;
		      G4int Z_new = Z;
		      
		      if (ptype == 1) {
			QPP_new--;
			Z_new--;
		      };
		      
		      if (ptype == 2) QNP_new--;
		      
		      if (verboseLevel > 3) {
			G4cout << " nucleus   px " << PEX.px() << " py " << PEX.py()
			       << " pz " << PEX.pz() << " E " << PEX.e() << G4endl
			       << " evaporate px " << mom.px() << " py " << mom.py()
			       << " pz " << mom.pz() << " E " << mom.e() << G4endl;
		      }
		
		      // New excitation energy depends on residual nuclear state
		      G4double mass_new = G4InuclNuclei::getNucleiMass(A_new, Z_new);
		      
		      EEXS_new = ((PEX-mom).m() - mass_new)*GeV;
		      if (EEXS_new < 0.) continue;	// Sanity check for new nucleus
		      
		      if (verboseLevel > 3)
			G4cout << " EEXS_new " << EEXS_new << G4endl;
		      
		      PEX -= mom;
		      EEXS = EEXS_new;
		      
		      A = A_new;
		      Z = Z_new;
		      
		      NEX--;
		      QEX--;
		      QP--;
		      QPP = QPP_new;
		      QNP = QNP_new;
		      
		      particle.setMomentum(mom);
		      output.addOutgoingParticle(particle);
		      ppout += mom;
		      if (verboseLevel > 3) {
			G4cout << particle << G4endl
			       << " ppout px " << ppout.px() << " py " << ppout.py()
			       << " pz " << ppout.pz() << " E " << ppout.e() << G4endl;
		      }

		      bad = false;
		    }		// while (itry1<itry_max && icase>0
		    
		    if (itry1 == itry_max) icase = 0;
		  }	// if (E < 0.) [else]
		}	// if (try_again && icase != 0)
	      }		// if (try_again)
	    }		// if (icase > 0)
	    
	    if (icase == 0 && try_again) { // N -> N + 2 
	      if (verboseLevel > 3) G4cout << " adding excitons" << G4endl;

	      G4double TNN = 1.6 * EFN + ESP;
	      G4double TNP = 1.6 * EFP + ESP;
	      G4double XNUN = 1.0 / (1.6 + ESP / EFN);
	      G4double XNUP = 1.0 / (1.6 + ESP / EFP);
	      G4double SNN1 = csNN(TNP) * XNUP;
	      G4double SNN2 = csNN(TNN) * XNUN;
	      G4double SPN1 = csPN(TNP) * XNUP;
	      G4double SPN2 = csPN(TNN) * XNUN;
	      G4double PP = (QPP * SNN1 + QNP * SPN1) * ZR;
	      G4double PN = (QPP * SPN2 + QNP * SNN2) * (AR - ZR);
	      G4double PW = PP + PN;
	      NEX += 2;
	      QEX += 2; 
	      QP++;
	      QH++;
	      AR--;
	      
	      if (AR > 1) {
		G4double SL = PW * inuclRndm();
		
		if (SL > PP) {
		  QNP++;
		  QNH++;
		} else {
		  QPP++;
		  QPH++;
		  ZR--;
		  if (ZR < 2) try_again = false;
		}  
	      } else try_again = false;
	    }	// if (icase==0 && try_again)
	  }	// if (try_again)
	} else try_again = false;	// if (EMN > eexs_cut)
      } else try_again = false;		// if (QEX < sqrg(2*EG)
    } else try_again = false;		// if (A > a_cut ...
  }		// while (try_again)
  
  // everything finished, set output fragment

  if (output.numberOfOutgoingParticles() == 0) {
    output.addRecoilFragment(target);
  } else {
    G4LorentzVector pnuc = pin - ppout;
    output.addRecoilFragment(makeFragment(pnuc, A, Z, EEXS));
    
    if (verboseLevel>3) 
      G4cout << " remaining nucleus\n" << output.getRecoilFragment() << G4endl;
  }

  validateOutput(target, output);	// Check energy conservation, etc.
  return;
}

G4double G4NonEquilibriumEvaporator::getMatrixElement(G4int a) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4NonEquilibriumEvaporator::getMatrixElement" << G4endl;
  }

  G4double me;

  if (a > 150) me = 100.0;
  else if (a > 20) me = 140.0;
  else me = 70.0;
 
  return me;
}

G4double G4NonEquilibriumEvaporator::getE0(G4int ) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4NonEquilibriumEvaporator::getEO" << G4endl;
  }

  const G4double e0 = 200.0;

  return e0;   
}

G4double G4NonEquilibriumEvaporator::getParLev(G4int a, G4int ) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4NonEquilibriumEvaporator::getParLev" << G4endl;
  }

  //  const G4double par = 0.125;
  G4double pl = 0.125 * a;

  return pl; 
}
