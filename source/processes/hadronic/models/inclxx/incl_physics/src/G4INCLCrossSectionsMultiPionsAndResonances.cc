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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLCrossSectionsMultiPionsAndResonances.cc
 * \brief Multipion and mesonic Resonances cross sections
 *
 * \date 4th February 2014
 * \author Jean-Christophe David
 */

#include "G4INCLCrossSectionsMultiPionsAndResonances.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
// #include <cassert>

namespace G4INCL {
	
	template<G4int N>
    struct BystrickyEvaluator {
		static G4double eval(const G4double pLab, const G4double oneOverThreshold, HornerCoefficients<N> const &coeffs) {
			const G4double pMeV = pLab*1E3;
			const G4double ekin=std::sqrt(ParticleTable::effectiveNucleonMass2+pMeV*pMeV)-ParticleTable::effectiveNucleonMass;
			const G4double xrat=ekin*oneOverThreshold;
			const G4double x=std::log(xrat);
			return HornerEvaluator<N>::eval(x, coeffs) * x * std::exp(-0.5*x);
		}
    };
	
	const G4int CrossSectionsMultiPionsAndResonances::nMaxPiNN = 4;
	const G4int CrossSectionsMultiPionsAndResonances::nMaxPiPiN = 4;
	
	const G4double CrossSectionsMultiPionsAndResonances::s11pzOOT = 0.0035761542037692665889;
	const G4double CrossSectionsMultiPionsAndResonances::s01ppOOT = 0.003421025623481919853;
	const G4double CrossSectionsMultiPionsAndResonances::s01pzOOT = 0.0035739814152966403123;
	const G4double CrossSectionsMultiPionsAndResonances::s11pmOOT = 0.0034855350296270480281;
	const G4double CrossSectionsMultiPionsAndResonances::s12pmOOT = 0.0016672224074691565119;
	const G4double CrossSectionsMultiPionsAndResonances::s12ppOOT = 0.0016507643038726931312;
	const G4double CrossSectionsMultiPionsAndResonances::s12zzOOT = 0.0011111111111111111111;
	const G4double CrossSectionsMultiPionsAndResonances::s02pzOOT = 0.00125;
	const G4double CrossSectionsMultiPionsAndResonances::s02pmOOT = 0.0016661112962345883443;
	const G4double CrossSectionsMultiPionsAndResonances::s12mzOOT = 0.0017047391749062392793;
	
	CrossSectionsMultiPionsAndResonances::CrossSectionsMultiPionsAndResonances() :
    s11pzHC(-2.228000000000294018,8.7560000000005723725,-0.61000000000023239325,-5.4139999999999780324,3.3338333333333348023,-0.75835000000000022049,0.060623611111111114688),
    s01ppHC(2.0570000000126518344,-6.029000000012135826,36.768500000002462784,-45.275666666666553533,25.112666666666611953,-7.2174166666666639187,1.0478875000000000275,-0.060804365079365080846),
    s01pzHC(0.18030000000000441851,7.8700999999999953598,-4.0548999999999990425,0.555199999999999959),
    s11pmHC(0.20590000000000031866,3.3450999999999993936,-1.4401999999999997825,0.17076666666666664973),
    s12pmHC(-0.77235999999999901328,4.2626599999999991117,-1.9008899999999997323,0.30192266666666663379,-0.012270833333333331986),
    s12ppHC(-0.75724999999999975664,2.0934399999999998565,-0.3803099999999999814),
    s12zzHC(-0.89599999999996965072,7.882999999999978632,-7.1049999999999961928,1.884333333333333089),
    s02pzHC(-1.0579999999999967036,11.113999999999994089,-8.5259999999999990196,2.0051666666666666525),
    s02pmHC(2.4009000000012553286,-7.7680000000013376183,20.619000000000433505,-16.429666666666723928,5.2525708333333363472,-0.58969166666666670206),
    s12mzHC(-0.21858699999999976269,1.9148999999999999722,-0.31727500000000001065,-0.027695000000000000486)
	{
	}
	
	G4double CrossSectionsMultiPionsAndResonances::total(Particle const * const p1, Particle const * const p2) {
		G4double inelastic;
		if(p1->isNucleon() && p2->isNucleon()) {
			return CrossSectionsMultiPions::NNTot(p1, p2);
		} else if((p1->isNucleon() && p2->isDelta()) ||
				  (p1->isDelta() && p2->isNucleon())) {
			inelastic = CrossSectionsMultiPions::NDeltaToNN(p1, p2);
		} else if((p1->isNucleon() && p2->isPion()) ||
				  (p1->isPion() && p2->isNucleon())) {
			return CrossSectionsMultiPions::piNTot(p1,p2);
		} else if((p1->isNucleon() && p2->isEta()) ||
				  (p1->isEta() && p2->isNucleon())) {
			inelastic = etaNToPiN(p1,p2) + etaNToPiPiN(p1,p2);
		} else if((p1->isNucleon() && p2->isOmega()) ||
				  (p1->isOmega() && p2->isNucleon())) {
			inelastic = omegaNInelastic(p1,p2);
		} else if((p1->isNucleon() && p2->isEtaPrime()) ||
				  (p1->isEtaPrime() && p2->isNucleon())) {
			inelastic = etaPrimeNToPiN(p1,p2);
		} else {
			inelastic = 0.;
		}
		
		return inelastic + elastic(p1, p2);
	}	

	
	G4double CrossSectionsMultiPionsAndResonances::elastic(Particle const * const p1, Particle const * const p2) {
		if((p1->isNucleon()||p1->isDelta()) && (p2->isNucleon()||p2->isDelta())){
			return CrossSectionsMultiPions::elastic(p1, p2);
		}
		else if ((p1->isNucleon() && p2->isPion()) || (p2->isNucleon() && p1->isPion())){
			return CrossSectionsMultiPions::elastic(p1, p2);
		}
		else if ((p1->isNucleon() && p2->isEta()) || (p2->isNucleon() && p1->isEta())){
			return etaNElastic(p1, p2);
		}
		else if ((p1->isNucleon() && p2->isOmega()) || (p2->isNucleon() && p1->isOmega())){
			return omegaNElastic(p1, p2);
		}
		else {
			return 0.0;
		}
	}
	
	
	G4double CrossSectionsMultiPionsAndResonances::piNToxPiN(const G4int xpi, Particle const * const particle1, Particle const * const particle2) {
			//
			//     pion-Nucleon producing xpi pions cross sections (corrected due to eta and omega)
			//
// assert(xpi>1 && xpi<=nMaxPiPiN);
// assert((particle1->isNucleon() && particle2->isPion()) || (particle1->isPion() && particle2->isNucleon()));

			const G4double oldXS2Pi=CrossSectionsMultiPions::piNToxPiN(2,particle1, particle2);
			const G4double oldXS3Pi=CrossSectionsMultiPions::piNToxPiN(3,particle1, particle2);
			const G4double oldXS4Pi=CrossSectionsMultiPions::piNToxPiN(4,particle1, particle2);
			const G4double xsEta=piNToEtaN(particle1, particle2);
			const G4double xsOmega=piNToOmegaN(particle1, particle2);
			G4double newXS2Pi=0.;	
			G4double newXS3Pi=0.;	
			G4double newXS4Pi=0.;
			
			if (xpi == 2) {
				if (oldXS4Pi != 0.)
					newXS2Pi=oldXS2Pi;
				else if (oldXS3Pi != 0.) {
					newXS3Pi=oldXS3Pi-xsEta-xsOmega;
					if (newXS3Pi < 1.e-09)
						newXS2Pi=oldXS2Pi-(xsEta+xsOmega-oldXS3Pi);
					else
						newXS2Pi=oldXS2Pi;
				}
				else { 
					newXS2Pi=oldXS2Pi-xsEta-xsOmega;
					if (newXS2Pi < 1.e-09)
						newXS2Pi=0.;
				}
				return newXS2Pi;
			}											
			else if (xpi == 3) {
				if (oldXS4Pi != 0.) {
					newXS4Pi=oldXS4Pi-xsEta-xsOmega;
					if (newXS4Pi < 1.e-09)
						newXS3Pi=oldXS3Pi-(xsEta+xsOmega-oldXS4Pi);
					else
						newXS3Pi=oldXS3Pi;
				}
				else { 
					newXS3Pi=oldXS3Pi-xsEta-xsOmega;
					if (newXS3Pi < 1.e-09)
						newXS3Pi=0.;
				}
				return newXS3Pi;
			}									
			else if (xpi == 4) {
				newXS4Pi=oldXS4Pi-xsEta-xsOmega;
				if (newXS4Pi < 1.e-09)
					newXS4Pi=0.;
				return newXS4Pi;
			}			
			else // should never reach this point
				return 0.;
		}
	
	G4double CrossSectionsMultiPionsAndResonances::piNToEtaN(Particle const * const particle1, Particle const * const particle2) {
		//
		//     Pion-Nucleon producing Eta cross sections
		//
// assert((particle1->isNucleon() && particle2->isPion()) || (particle1->isPion() && particle2->isNucleon()));
		
		G4double sigma;
		sigma=piMinuspToEtaN(particle1,particle2);
		
		const G4int isoin = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
		
		if (isoin == -1) {
			if ((particle1->getType()) == Proton || (particle2->getType()) == Proton) return sigma;
			else return 0.5 * sigma;
		}
		else if (isoin == 1) {
			if ((particle1->getType()) == Neutron || (particle2->getType()) == Neutron) return sigma;
			else return 0.5 * sigma;
		}
		else return 0. ; // should never return 0. (?) // pi+ p and pi- n return 0.
		
//		return sigma;
	}
	
	G4double CrossSectionsMultiPionsAndResonances::piNToOmegaN(Particle const * const particle1, Particle const * const particle2) {
		//
		//     Pion-Nucleon producing Omega cross sections
		//
// assert((particle1->isNucleon() && particle2->isPion()) || (particle1->isPion() && particle2->isNucleon()));
		
		G4double sigma;
		sigma=piMinuspToOmegaN(particle1,particle2);
		
		const G4int isoin = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
		
  if (isoin == -1) {
			if ((particle1->getType()) == Proton || (particle2->getType()) == Proton) return sigma;
			else return 0.5 * sigma;
		}
		else if (isoin == 1) {
			if ((particle1->getType()) == Neutron || (particle2->getType()) == Neutron) return sigma;
			else return 0.5 * sigma;
		}
		else return 0. ; // should never return 0. (?) // pi+ p and pi- n return 0.
		
//		return sigma;
	}
	
#if defined(NDEBUG) || defined(INCLXX_IN_GEANT4_MODE)
    G4double CrossSectionsMultiPionsAndResonances::piNToEtaPrimeN(Particle const * const /*particle1*/, Particle const * const /*particle2*/) {
#else
	G4double CrossSectionsMultiPionsAndResonances::piNToEtaPrimeN(Particle const * const particle1, Particle const * const particle2) {
#endif
		//
		//     Pion-Nucleon producing EtaPrime cross sections
		//
// assert((particle1->isNucleon() && particle2->isPion()) || (particle1->isPion() && particle2->isNucleon()));
		
		return 0.;
	}
	
	G4double CrossSectionsMultiPionsAndResonances::etaNToPiN(Particle const * const particle1, Particle const * const particle2) {
 //
 //     Eta-Nucleon producing Pion cross sections
 //
// assert((particle1->isNucleon() && particle2->isEta()) || (particle1->isEta() && particle2->isNucleon()));
 
 const Particle *eta;
 const Particle *nucleon;
 
 if (particle1->isEta()) {
  eta = particle1;
  nucleon = particle2;
 }
 else {
  eta = particle2;
  nucleon = particle1;
 }
 
 const G4double pLab = KinematicsUtils::momentumInLab(eta, nucleon);
 G4double sigma=0.;		
 
 if (pLab <= 574.) 
  sigma= 1.511147E-13*std::pow(pLab,6)- 3.603636E-10*std::pow(pLab,5)+ 3.443487E-07*std::pow(pLab,4)- 1.681980E-04*std::pow(pLab,3)+ 4.437913E-02*std::pow(pLab,2)- 6.172108E+00*pLab+ 4.031449E+02;
 else if (pLab <= 850.) 
  sigma= -8.00018E-14*std::pow(pLab,6)+ 3.50041E-10*std::pow(pLab,5)- 6.33891E-07*std::pow(pLab,4)+ 6.07658E-04*std::pow(pLab,3)- 3.24936E-01*std::pow(pLab,2)+ 9.18098E+01*pLab- 1.06943E+04;
 else if (pLab <= 1300.)
  sigma= 6.56364E-09*std::pow(pLab,3)- 2.07653E-05*std::pow(pLab,2)+ 1.84148E-02*pLab- 1.70427E+00;	 
 else {
  G4double ECM=KinematicsUtils::totalEnergyInCM(eta, nucleon);
		G4double massPiZero=ParticleTable::getINCLMass(PiZero);
		G4double massPiMinus=ParticleTable::getINCLMass(PiMinus);
		G4double massProton=ParticleTable::getINCLMass(Proton);
		G4double masseta;
		G4double massnucleon;
		G4double pCM_eta;
  masseta=eta->getMass();
  massnucleon=nucleon->getMass();
  pCM_eta=KinematicsUtils::momentumInCM(ECM, masseta, massnucleon);		 
  G4double pCM_PiZero=KinematicsUtils::momentumInCM(ECM, massPiZero, massProton);
  G4double pCM_PiMinus=KinematicsUtils::momentumInCM(ECM, massPiMinus, massProton); // = pCM_PiPlus (because massPiMinus = massPiPlus) 
  sigma = (piMinuspToEtaN(ECM)/2.) * std::pow((pCM_PiZero/pCM_eta), 2) + piMinuspToEtaN(ECM) * std::pow((pCM_PiMinus/pCM_eta), 2);			
 } 
 if (sigma < 0.) sigma=0.;
 return sigma;
}

	G4double CrossSectionsMultiPionsAndResonances::etaNToPiPiN(Particle const * const particle1, Particle const * const particle2) {
			//
			//     Eta-Nucleon producing Two Pions cross sections
			//
// assert((particle1->isNucleon() && particle2->isEta()) || (particle1->isEta() && particle2->isNucleon()));
						
			G4double sigma=0.;		
			
			const Particle *eta;
			const Particle *nucleon;

			if (particle1->isEta()) {
				eta = particle1;
				nucleon = particle2;
			}
		 else {
				eta = particle2;
				nucleon = particle1;
			}
				
		 const G4double pLab = KinematicsUtils::momentumInLab(eta, nucleon);
			
   if (pLab < 450.) 
				sigma = 2.01854221E-13*std::pow(pLab,6) - 3.49750459E-10*std::pow(pLab,5) + 2.46011585E-07*std::pow(pLab,4) - 9.01422901E-05*std::pow(pLab,3) + 1.83382964E-02*std::pow(pLab,2) - 2.03113098E+00*pLab + 1.10358550E+02;				
			else if (pLab < 600.) 
				sigma = 2.01854221E-13*std::pow(450.,6) - 3.49750459E-10*std::pow(450.,5) + 2.46011585E-07*std::pow(450.,4) - 9.01422901E-05*std::pow(450.,3) + 1.83382964E-02*std::pow(450.,2) - 2.03113098E+00*450. + 1.10358550E+02;				
			else if (pLab <= 1300.) 
				sigma = -6.32793049e-16*std::pow(pLab,6) + 3.95985900e-12*std::pow(pLab,5) - 1.01727714e-8*std::pow(pLab,4) + 
				         1.37055547e-05*std::pow(pLab,3) - 1.01830486e-02*std::pow(pLab,2) + 3.93492126*pLab - 609.447145;
   else 
				sigma = etaNToPiN(particle1,particle2);
			
			if (sigma < 0.) sigma = 0.;
			return sigma; // Parameterization from the ANL-Osaka DCC model [PRC88(2013)035209] - eta p --> "pi+pi0 n" + "pi0 pi0 p" total XS
		}
		
		
	G4double CrossSectionsMultiPionsAndResonances::etaNElastic(Particle const * const particle1, Particle const * const particle2) {
			//
			//     Eta-Nucleon elastic cross sections
			//
// assert((particle1->isNucleon() && particle2->isEta()) || (particle1->isEta() && particle2->isNucleon()));
			
			G4double sigma=0.;		
			
			const Particle *eta;
			const Particle *nucleon;
			
			if (particle1->isEta()) {
				eta = particle1;
				nucleon = particle2;
			}
		 else {
				eta = particle2;
				nucleon = particle1;
			}
			
		 const G4double pLab = KinematicsUtils::momentumInLab(eta, nucleon);
			
   if (pLab < 700.) 
				sigma = 3.6838e-15*std::pow(pLab,6) - 9.7815e-12*std::pow(pLab,5) + 9.7914e-9*std::pow(pLab,4) - 
				        4.3222e-06*std::pow(pLab,3) + 7.9188e-04*std::pow(pLab,2) - 1.8379e-01*pLab + 84.965;			
			else if (pLab < 1400.) 
				sigma = 3.562630e-16*std::pow(pLab,6) - 2.384766e-12*std::pow(pLab,5) + 6.601312e-9*std::pow(pLab,4) - 
				        9.667078e-06*std::pow(pLab,3) + 7.894845e-03*std::pow(pLab,2) - 3.409200*pLab + 609.8501;
			else if (pLab < 2025.)
				sigma = -1.041950E-03*pLab + 2.110529E+00;
			else	
			 sigma=0.;
				
			if (sigma < 0.) sigma = 0.;
			return sigma; // Parameterization from the ANL-Osaka DCC model [PRC88(2013)035209]
		}

	G4double CrossSectionsMultiPionsAndResonances::omegaNInelastic(Particle const * const particle1, Particle const * const particle2) {
		//
		//     Omega-Nucleon inelastic cross sections
		//
// assert((particle1->isNucleon() && particle2->isOmega()) || (particle1->isOmega() && particle2->isNucleon()));
		
		G4double sigma=0.;		
		
		const Particle *omega;
		const Particle *nucleon;
		
		if (particle1->isOmega()) {
			omega = particle1;
			nucleon = particle2;
			}
		else {
			omega = particle2;
			nucleon = particle1;
		}
		
		const G4double pLab = KinematicsUtils::momentumInLab(omega, nucleon)/1000.; // GeV/c
		
		sigma = 20. + 4.0/pLab;	// Eq.(24) in G.I. Lykasov et al., EPJA 6, 71-81 (1999)
		
		return sigma;
	}

  
	G4double CrossSectionsMultiPionsAndResonances::omegaNElastic(Particle const * const particle1, Particle const * const particle2) {
		//
		//     Omega-Nucleon elastic cross sections
		//
// assert((particle1->isNucleon() && particle2->isOmega()) || (particle1->isOmega() && particle2->isNucleon()));
		
		G4double sigma=0.;		
		
		const Particle *omega;
		const Particle *nucleon;
		
		if (particle1->isOmega()) {
			omega = particle1;
			nucleon = particle2;
		}
		else {
			omega = particle2;
			nucleon = particle1;
		}
		
		const G4double pLab = KinematicsUtils::momentumInLab(omega, nucleon)/1000.; // GeV/c
		
		sigma = 5.4 + 10.*std::exp(-0.6*pLab);	// Eq.(21) in G.I. Lykasov et al., EPJA 6, 71-81 (1999)
		
		return sigma;
	}
  
		
	G4double CrossSectionsMultiPionsAndResonances::omegaNToPiN(Particle const * const particle1, Particle const * const particle2) {
		//
		//     Omega-Nucleon producing Pion cross sections
		//
// assert((particle1->isNucleon() && particle2->isOmega()) || (particle1->isOmega() && particle2->isNucleon()));
		
		G4double ECM=KinematicsUtils::totalEnergyInCM(particle1, particle2);
		
		G4double massPiZero=ParticleTable::getINCLMass(PiZero);
		G4double massPiMinus=ParticleTable::getINCLMass(PiMinus);
		G4double massProton=ParticleTable::getINCLMass(Proton);
		
		G4double massomega;
		G4double massnucleon;
		G4double pCM_omega;
		G4double pLab_omega;
		
		G4double sigma=0.;		
		
		if (particle1->isOmega()) {
			massomega=particle1->getMass();
			massnucleon=particle2->getMass();
		}
		else {
			massomega=particle2->getMass();
			massnucleon=particle1->getMass();
		}
		pCM_omega=KinematicsUtils::momentumInCM(ECM, massomega, massnucleon);		
		pLab_omega=KinematicsUtils::momentumInLab(ECM*ECM, massomega, massnucleon);		
		
		G4double pCM_PiZero=KinematicsUtils::momentumInCM(ECM, massPiZero, massProton);
		G4double pCM_PiMinus=KinematicsUtils::momentumInCM(ECM, massPiMinus, massProton); // = pCM_PiPlus (because massPiMinus = massPiPlus)
		
		sigma = (piMinuspToOmegaN(ECM)/2.) * std::pow((pCM_PiZero/pCM_omega), 2) 
		+ piMinuspToOmegaN(ECM) * std::pow((pCM_PiMinus/pCM_omega), 2);			
  
  if (sigma > omegaNInelastic(particle1, particle2) || (pLab_omega < 200.)) {
//  if (sigma > omegaNInelastic(particle1, particle2)) {
   sigma = omegaNInelastic(particle1, particle2);
  }  
  
		return sigma;
	}
	
  
	G4double CrossSectionsMultiPionsAndResonances::omegaNToPiPiN(Particle const * const particle1, Particle const * const particle2) {
		//
		//     Omega-Nucleon producing 2 PionS cross sections
		//
// assert((particle1->isNucleon() && particle2->isOmega()) || (particle1->isOmega() && particle2->isNucleon()));
		
		G4double sigma=0.;		
		
		sigma = omegaNInelastic(particle1,particle2) - omegaNToPiN(particle1,particle2) ;
		
		return sigma;
	}
  
	
#if defined(NDEBUG) || defined(INCLXX_IN_GEANT4_MODE)
  G4double CrossSectionsMultiPionsAndResonances::etaPrimeNToPiN(Particle const * const /*particle1*/, Particle const * const /*particle2*/) {
#else
  G4double CrossSectionsMultiPionsAndResonances::etaPrimeNToPiN(Particle const * const particle1, Particle const * const particle2) {
#endif
		//
		//     EtaPrime-Nucleon producing Pion cross sections
		//
// assert((particle1->isNucleon() && particle2->isEtaPrime()) || (particle1->isEtaPrime() && particle2->isNucleon()));
		
		return 0.;
	}	
	
	G4double CrossSectionsMultiPionsAndResonances::piMinuspToEtaN(Particle const * const particle1, Particle const * const particle2) {
		//
		//     Pion-Nucleon producing Eta cross sections
		//
// assert((particle1->isNucleon() && particle2->isPion()) || (particle1->isPion() && particle2->isNucleon()));
		
		G4double masspion;
		G4double massnucleon;
		if (particle1->isPion()) {
		    masspion=particle1->getMass();
		    massnucleon=particle2->getMass();
		}
		else {
		    masspion=particle2->getMass();
		    massnucleon=particle1->getMass();
		}
		
		G4double ECM=KinematicsUtils::totalEnergyInCM(particle1, particle2);
		G4double plab=KinematicsUtils::momentumInLab(ECM*ECM, masspion, massnucleon)/1000.; // GeV/c
		
		G4double sigma;
		
// new parameterization (JCD) - end of february 2016
		if (ECM < 1486.5) sigma=0.;
		else
		{
			if (ECM < 1535.)
			{
				sigma = -0.0000003689197974814*std::pow(ECM,4) + 0.002260193900097*std::pow(ECM,3) - 5.193105877187*std::pow(ECM,2) + 5303.505273919*ECM - 2031265.900648;
			}
			else if (ECM < 1670.)
			{
				sigma = -0.0000000337986446*std::pow(ECM,4) + 0.000218279989*std::pow(ECM,3) - 0.528276144*std::pow(ECM,2) + 567.828367*ECM - 228709.42;
			}
			else if (ECM < 1714.)
			{
				sigma =  0.000003737765*std::pow(ECM,2) - 0.005664062*ECM;
			}
			else sigma=1.47*std::pow(plab, -1.68);
		}		
//		

		return sigma;
	}
	
	G4double CrossSectionsMultiPionsAndResonances::piMinuspToEtaN(const G4double ECM) {
		//
		//     Pion-Nucleon producing Eta cross sections
		//
		
  const G4double masspion = ParticleTable::getRealMass(PiMinus);		
  const G4double massnucleon = ParticleTable::getRealMass(Proton);		
  
  G4double plab=KinematicsUtils::momentumInLab(ECM*ECM, masspion, massnucleon)/1000.;  // GeV/c
		
		G4double sigma;
		
// new parameterization (JCD) - end of february 2016
		if (ECM < 1486.5) sigma=0.;
		else
		{
			if (ECM < 1535.)
			{
				sigma = -0.0000003689197974814*std::pow(ECM,4) + 0.002260193900097*std::pow(ECM,3) - 5.193105877187*std::pow(ECM,2) + 5303.505273919*ECM - 2031265.900648;
			}
			else if (ECM < 1670.)
			{
				sigma = -0.0000000337986446*std::pow(ECM,4) + 0.000218279989*std::pow(ECM,3) - 0.528276144*std::pow(ECM,2) + 567.828367*ECM - 228709.42;
			}
			else if (ECM < 1714.)
			{
				sigma =  0.000003737765*std::pow(ECM,2) - 0.005664062*ECM;
			}
			else sigma=1.47*std::pow(plab, -1.68);
		}		
//		

		return sigma;
	}

	G4double CrossSectionsMultiPionsAndResonances::piMinuspToOmegaN(Particle const * const particle1, Particle const * const particle2) {
		//
		//     Pion-Nucleon producing Omega cross sections
		//
// assert((particle1->isNucleon() && particle2->isPion()) || (particle1->isPion() && particle2->isNucleon()));
//jcd to be removed
//  return 0.;
//jcd    
		
//		G4double param=1.095 ; // Deneye (Thesis)
		G4double param=1.0903 ; // JCD (threshold taken into account)
  
  G4double masspion;
		G4double massnucleon;
		if (particle1->isPion()) {
		    masspion=particle1->getMass();
		    massnucleon=particle2->getMass();
		}
		else {
		    masspion=particle2->getMass();
		    massnucleon=particle1->getMass();
		}
		G4double ECM=KinematicsUtils::totalEnergyInCM(particle1, particle2);
		G4double plab=KinematicsUtils::momentumInLab(ECM*ECM, masspion, massnucleon)/1000.;  // GeV/c
		
		G4double sigma;
		if (plab < param) sigma=0.;
		else sigma=13.76*(plab-param)/(std::pow(plab, 3.33) - 1.07); // Phys. Rev. C 41, 1701â1718 (1990)
  
		return sigma;
}
			G4double CrossSectionsMultiPionsAndResonances::piMinuspToOmegaN(const G4double ECM) {
				//
				//     Pion-Nucleon producing Omega cross sections
				//
//jcd to be removed
//    return 0.;
//jcd    
				
//		G4double param=1.095 ; // Deneye (Thesis)
    G4double param=1.0903 ; // JCD (threshold taken into account)

				const G4double masspion = ParticleTable::getRealMass(PiMinus);		
				const G4double massnucleon = ParticleTable::getRealMass(Proton);		

				G4double plab=KinematicsUtils::momentumInLab(ECM*ECM, masspion, massnucleon)/1000.;  // GeV/c
				
				G4double sigma;
				if (plab < param) sigma=0.;
				else sigma=13.76*(plab-param)/(std::pow(plab, 3.33)-1.07);
				
				return sigma;
}

    G4double CrossSectionsMultiPionsAndResonances::NNToNNEtaIso(const G4double ener, const G4int iso) {
		
		const G4double Ecm=0.001*ener;
		G4double sNNEta; // pp->pp+eta(+X)
		G4double sNNEta1; // np->np+eta(+X)
		G4double sNNEta2; // np->d+eta (d will be considered as np - How far is this right?)
		G4double x=Ecm*Ecm/5.88;
                
//jcd
  if (Ecm >= 3.05) {
     sNNEta = 2.5*std::pow((x-1.),1.47)*std::pow(x,-1.25)*1000.;
  }      
  else if(Ecm >= 2.6) {
     sNNEta = -327.29*Ecm*Ecm*Ecm + 2870.*Ecm*Ecm - 7229.3*Ecm + 5273.3;
     if (sNNEta <= NNToNNEtaExcluIso(ener, 2)*1000.) sNNEta = NNToNNEtaExcluIso(ener, 2)*1000.;
  }
  else {
     sNNEta = NNToNNEtaExcluIso(ener, 2)*1000.;
  }
//jcd
  if (sNNEta < 1.e-9) sNNEta = 0.;

  if (iso != 0) {
     return sNNEta/1000.; // parameterization in microbarn (not millibarn)!
  }

  if(Ecm >= 6.25) {
     sNNEta1 = sNNEta;
  }
  else if (Ecm >= 2.6) {
     sNNEta1 = sNNEta*std::exp(-(-5.53151576/Ecm+0.8850425));
  }
  else if (Ecm >= 2.525) { // = exclusive pn
      sNNEta1 = -4433.586*Ecm*Ecm*Ecm*Ecm + 56581.54*Ecm*Ecm*Ecm - 270212.6*Ecm*Ecm + 571650.6*Ecm - 451091.6;
  }
  else { // = exclusive pn
      sNNEta1 = 17570.217219*Ecm*Ecm - 84910.985402*Ecm + 102585.55847;
  }

  sNNEta2 = -10220.89518466*Ecm*Ecm+51227.30841724*Ecm-64097.96025731;
  if (sNNEta2 < 0.) sNNEta2=0.;

 	sNNEta = 2*(sNNEta1+sNNEta2)-sNNEta;

  G4double Mn=ParticleTable::getRealMass(Neutron)/1000.;
  G4double Mp=ParticleTable::getRealMass(Proton)/1000.;
  G4double Meta=ParticleTable::getRealMass(Eta)/1000.;
		if (sNNEta < 1.e-9 || Ecm < Mn+Mp+Meta) sNNEta = 0.;
		
		return sNNEta/1000.; // parameterization in microbarn (not millibarn)!
	}

	
  G4double CrossSectionsMultiPionsAndResonances::NNToNNEta(Particle const * const particle1, Particle const * const particle2) {
		
  const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2);
  const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
		
		if (iso != 0) {
   return NNToNNEtaIso(ener, iso);
		}
		else {
			return 0.5*(NNToNNEtaIso(ener, 0)+NNToNNEtaIso(ener, 2));
		}
    }
	
  G4double CrossSectionsMultiPionsAndResonances::NNToNNEtaExcluIso(const G4double ener, const G4int iso) {
				
  const G4double Ecm=0.001*ener;
  G4double sNNEta; // pp->pp+eta
  G4double sNNEta1; // np->np+eta
  G4double sNNEta2; // np->d+eta (d wil be considered as np - How far is this right?)
				
  if(Ecm>=3.875) { // By hand (JCD)
     sNNEta = -13.008*Ecm*Ecm + 84.531*Ecm + 36.234;
  }
  else if(Ecm>=2.725) { // By hand (JCD)
     sNNEta = -913.2809*std::pow(Ecm,5) + 15564.27*std::pow(Ecm,4) - 105054.9*std::pow(Ecm,3) + 351294.2*std::pow(Ecm,2) - 582413.9*Ecm + 383474.7;
  }
  else if(Ecm>=2.575) { // By hand (JCD)
     sNNEta = -2640.3*Ecm*Ecm + 14692*Ecm - 20225;
  }
  else {
     sNNEta = -147043.497285*std::pow(Ecm,4) + 1487222.5438123*std::pow(Ecm,3) - 5634399.900744*std::pow(Ecm,2) + 9477290.199378*Ecm - 5972174.353438;
  }
				
   G4double Mn=ParticleTable::getRealMass(Neutron)/1000.;
   G4double Mp=ParticleTable::getRealMass(Proton)/1000.;
   G4double Meta=ParticleTable::getRealMass(Eta)/1000.;
   G4double Thr0=0.;
   if (iso > 0) {
    Thr0=2.*Mp+Meta;
   }
   else if (iso < 0) {
    Thr0=2.*Mn+Meta;
   }
   else {
    Thr0=Mn+Mp+Meta;
   }

   if (sNNEta < 1.e-9 || Ecm < Thr0) sNNEta = 0.;  // Thr0: Ecm threshold
				
   if (iso != 0) {
      return sNNEta/1000.; // parameterization in microbarn (not millibarn)!
   }
				
   if(Ecm>=3.9) {
      sNNEta1 = sNNEta;
   }
   else if (Ecm >= 3.5) {
      sNNEta1 = -1916.2*Ecm*Ecm*Ecm + 21556.0*Ecm*Ecm - 80828.0*Ecm + 101200.0;
   }
   else if (Ecm >= 2.525) {
      sNNEta1 = -4433.586*Ecm*Ecm*Ecm*Ecm + 56581.54*Ecm*Ecm*Ecm - 270212.6*Ecm*Ecm + 571650.6*Ecm - 451091.6;
   }
   else {
      sNNEta1 = 17570.217219*Ecm*Ecm - 84910.985402*Ecm + 102585.55847;
  }
    
   sNNEta2 = -10220.89518466*Ecm*Ecm+51227.30841724*Ecm-64097.96025731;
   if (sNNEta2 < 0.) sNNEta2=0.;
    
   sNNEta = 2*(sNNEta1+sNNEta2)-sNNEta;
   if (sNNEta < 1.e-9 || Ecm < Thr0) sNNEta = 0.;  // Thr0: Ecm threshold
				
  return sNNEta/1000.; // parameterization in microbarn (not millibarn)!

}
			
  G4double CrossSectionsMultiPionsAndResonances::NNToNNEtaExclu(Particle const * const particle1, Particle const * const particle2) {
				
   const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2);
   const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
				
   if (iso != 0) {
      return NNToNNEtaExcluIso(ener, iso);
   }
   else {
      return 0.5*(NNToNNEtaExcluIso(ener, 0)+NNToNNEtaExcluIso(ener, 2));
   }
}

   
   G4double CrossSectionsMultiPionsAndResonances::NNToNNOmegaIso(const G4double ener, const G4int iso) {
    
    const G4double Ecm=0.001*ener;
    G4double sNNOmega; // pp->pp+eta(+X)
    G4double sNNOmega1; // np->np+eta(+X)
    G4double x=Ecm*Ecm/7.06;
    
    if(Ecm>4.0) {
     sNNOmega = 2.5*std::pow(x-1, 1.47)*std::pow(x, -1.11) ;
    }
    else if(Ecm>2.802) { // 2802 MeV -> threshold to open inclusive (based on multipion threshold and omega mass) 
     sNNOmega = (568.5254*Ecm*Ecm - 2694.045*Ecm + 3106.247)/1000.;
     if (sNNOmega <= NNToNNOmegaExcluIso(ener, 2)) sNNOmega = NNToNNOmegaExcluIso(ener, 2);
    }
    else {
     sNNOmega = NNToNNOmegaExcluIso(ener, 2);
    }
    
    if (sNNOmega < 1.e-9) sNNOmega = 0.;
    
    if (iso != 0) {
     return sNNOmega; 
    }
    
    sNNOmega1 = 3.*sNNOmega; // 3.0: ratio pn/pp (5 from theory; 2 from experiments)
        
    sNNOmega = 2.*sNNOmega1-sNNOmega;
    
    if (sNNOmega < 1.e-9) sNNOmega = 0.;
    
    return sNNOmega; 
   }
   
   
   G4double CrossSectionsMultiPionsAndResonances::NNToNNOmega(Particle const * const particle1, Particle const * const particle2) {
    
    const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2);
    const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
//jcd to be removed
//    return 0.;
//jcd    
    if (iso != 0) {
     return NNToNNOmegaIso(ener, iso);
    }
    else {
     return 0.5*(NNToNNOmegaIso(ener, 0)+NNToNNOmegaIso(ener, 2));
    }
   }
   
   G4double CrossSectionsMultiPionsAndResonances::NNToNNOmegaExcluIso(const G4double ener, const G4int iso) {
				
				const G4double Ecm=0.001*ener;
				G4double sNNOmega; // pp->pp+eta
				G4double sNNOmega1; // np->np+eta
    G4double sthroot=std::sqrt(7.06);
				
				if(Ecm>=3.0744) { // By hand (JCD)
					sNNOmega = 330.*(Ecm-sthroot)/(1.05+std::pow((Ecm-sthroot),2));
				}
				else if(Ecm>=2.65854){
					sNNOmega = -1208.09757*std::pow(Ecm,3) + 10773.3322*std::pow(Ecm,2) - 31661.0223*Ecm + 30728.7241 ;
				}
				else {
					sNNOmega = 0. ;
				}
				
    G4double Mn=ParticleTable::getRealMass(Neutron)/1000.;
    G4double Mp=ParticleTable::getRealMass(Proton)/1000.;
    G4double Momega=ParticleTable::getRealMass(Omega)/1000.;
    G4double Thr0=0.;
    if (iso > 0) {
     Thr0=2.*Mp+Momega;
    }
    else if (iso < 0) {
     Thr0=2.*Mn+Momega;
    }
    else {
     Thr0=Mn+Mp+Momega;
    }
    
    if (sNNOmega < 1.e-9 || Ecm < Thr0) sNNOmega = 0.;  // Thr0: Ecm threshold
				
				if (iso != 0) {
					return sNNOmega/1000.; // parameterization in microbarn (not millibarn)!
				}
				
    sNNOmega1 = 3*sNNOmega; // 3.0: ratio pn/pp
        
				sNNOmega = 2*sNNOmega1-sNNOmega;
				if (sNNOmega < 1.e-9 || Ecm < Thr0) sNNOmega = 0.;
				
				return sNNOmega/1000.; // parameterization in microbarn (not millibarn)!
			}
			
   G4double CrossSectionsMultiPionsAndResonances::NNToNNOmegaExclu(Particle const * const particle1, Particle const * const particle2) {
//jcd to be removed
//    return 0.;
//jcd    
				
				const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2);
				const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
				
				if (iso != 0) {
					return NNToNNOmegaExcluIso(ener, iso);
				}
				else {
					return 0.5*(NNToNNOmegaExcluIso(ener, 0)+NNToNNOmegaExcluIso(ener, 2));
				}
			}
   
   
		G4double CrossSectionsMultiPionsAndResonances::NNToxPiNN(const G4int xpi, Particle const * const particle1, Particle const * const particle2) {
		//
		//     Nucleon-Nucleon producing xpi pions cross sections
		//
// assert(xpi>0 && xpi<=nMaxPiNN);
// assert(particle1->isNucleon() && particle2->isNucleon());
		
		G4double oldXS1Pi=CrossSectionsMultiPions::NNToxPiNN(1,particle1, particle2);
		G4double oldXS2Pi=CrossSectionsMultiPions::NNToxPiNN(2,particle1, particle2);
		G4double oldXS3Pi=CrossSectionsMultiPions::NNToxPiNN(3,particle1, particle2);
		G4double oldXS4Pi=CrossSectionsMultiPions::NNToxPiNN(4,particle1, particle2);
		G4double xsEtaOmega=NNToNNEta(particle1, particle2)+NNToNNOmega(particle1, particle2);
		G4double newXS1Pi=0.;	
		G4double newXS2Pi=0.;	
		G4double newXS3Pi=0.;	
		G4double newXS4Pi=0.;	
			
		if (xpi == 1) {
			if (oldXS4Pi != 0. || oldXS3Pi != 0.)
				newXS1Pi=oldXS1Pi;
			else if (oldXS2Pi != 0.) {
				newXS2Pi=oldXS2Pi-xsEtaOmega;
			 if (newXS2Pi < 0.)
				 newXS1Pi=oldXS1Pi-(xsEtaOmega-oldXS2Pi);
				else
					newXS1Pi=oldXS1Pi;
			}
			else 
				newXS1Pi=oldXS1Pi-xsEtaOmega;
		 return newXS1Pi;
		}
		else if (xpi == 2) {
			if (oldXS4Pi != 0.)
				newXS2Pi=oldXS2Pi;
			else if (oldXS3Pi != 0.) {
				newXS3Pi=oldXS3Pi-xsEtaOmega;
				if (newXS3Pi < 0.)
					newXS2Pi=oldXS2Pi-(xsEtaOmega-oldXS3Pi);
				else
					newXS2Pi=oldXS2Pi;
			}
			else { 
				newXS2Pi=oldXS2Pi-xsEtaOmega;
				if (newXS2Pi < 0.)
					newXS2Pi=0.;
			}
			return newXS2Pi;
		}
		else if (xpi == 3) {
			if (oldXS4Pi != 0.) {
				newXS4Pi=oldXS4Pi-xsEtaOmega;
			 if (newXS4Pi < 0.)
			 	newXS3Pi=oldXS3Pi-(xsEtaOmega-oldXS4Pi);
			 else
			 	newXS3Pi=oldXS3Pi;
			}
			else { 
				newXS3Pi=oldXS3Pi-xsEtaOmega;				
				if (newXS3Pi < 0.)
					newXS3Pi=0.;
			}
			return newXS3Pi;
		}
		else if (xpi == 4) {
			newXS4Pi=oldXS4Pi-xsEtaOmega;
			if (newXS4Pi < 0.)
				newXS4Pi=0.;
		return newXS4Pi;
		}

		else // should never reach this point
			return 0.;
	}


			G4double CrossSectionsMultiPionsAndResonances::NNToNNEtaOnePi(Particle const * const particle1, Particle const * const particle2) {
				// Cross section for nucleon-nucleon producing one eta and one pion
				
				const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
				if (iso!=0)
					return 0.;
				
				const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 581.437; // 581.437 MeV translation to open pion production in NNEta (= 2705.55 - 2018.563; 4074595.287720512986=2018.563*2018.563)
				if (ener < 2018.563) return 0.;
				
				const G4double xsiso2=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
				const G4double xsiso0=CrossSectionsMultiPions::NNInelasticIso(ener, 0);

				return 0.25*(CrossSectionsMultiPions::NNOnePiOrDelta(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNOnePiOrDelta(ener, 2, xsiso2));
			}
			
		 	G4double CrossSectionsMultiPionsAndResonances::NNToNNEtaOnePiOrDelta(Particle const * const particle1, Particle const * const particle2) {
				const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 581.437; // 581.437 MeV translation to open pion production in NNEta
				if (ener < 2018.563) return 0.;
				const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
				
				const G4double xsiso2=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
				if (iso != 0)
					return CrossSectionsMultiPions::NNOnePiOrDelta(ener, iso, xsiso2);
				else {
					const G4double xsiso0=CrossSectionsMultiPions::NNInelasticIso(ener, 0);
					return 0.5*(CrossSectionsMultiPions::NNOnePiOrDelta(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNOnePiOrDelta(ener, 2, xsiso2));
				}
			}
			
			G4double CrossSectionsMultiPionsAndResonances::NNToNNEtaTwoPi(Particle const * const particle1, Particle const * const particle2) {
				//
				//     Nucleon-Nucleon producing one eta and two pions 
				//
				const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 581.437; // 581.437 MeV translation to open pion production in NNEta
				if (ener < 2018.563) return 0.;
				const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
				
				
				const G4double xsiso2=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
				if (iso != 0) {
					return CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2);
				}
				else {
					const G4double xsiso0=CrossSectionsMultiPions::NNInelasticIso(ener, 0);
					return 0.5*(CrossSectionsMultiPions::NNTwoPi(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2));
				}
			}
			
			G4double CrossSectionsMultiPionsAndResonances::NNToNNEtaThreePi(Particle const * const particle1, Particle const * const particle2) {
				//
				//     Nucleon-Nucleon producing one eta and three pions
				//
				
				const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 581.437; // 581.437 MeV translation to open pion production in NNEta
				if (ener < 2018.563) return 0.;
				const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
				
				
				const G4double xsiso2=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
				const G4double xs1pi2=CrossSectionsMultiPions::NNOnePiOrDelta(ener, 2, xsiso2);
				const G4double xs2pi2=CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2);
				if (iso != 0)
					return CrossSectionsMultiPions::NNThreePi(ener, 2, xsiso2, xs1pi2, xs2pi2);
				else {
					const G4double xsiso0=CrossSectionsMultiPions::NNInelasticIso(ener, 0);
					const G4double xs1pi0=CrossSectionsMultiPions::NNOnePiOrDelta(ener, 0, xsiso0);
					const G4double xs2pi0=CrossSectionsMultiPions::NNTwoPi(ener, 0, xsiso0);
					return 0.5*(CrossSectionsMultiPions::NNThreePi(ener, 0, xsiso0, xs1pi0, xs2pi0)+ CrossSectionsMultiPions::NNThreePi(ener, 2, xsiso2, xs1pi2, xs2pi2));
				}
			}
			
			G4double CrossSectionsMultiPionsAndResonances::NNToNNEtaFourPi(Particle const * const particle1, Particle const * const particle2) {
				//
				//     Nucleon-Nucleon producing one eta and four pions
				//
				
				const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 581.437; // 581.437 MeV translation to open pion production in NNEta
				if (ener < 2018.563) return 0.;
				const G4double s = ener*ener;
				const G4int i = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
    G4double xsinelas;
				if (i!=0) 
     xsinelas=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
    else 
					xsinelas=0.5*(CrossSectionsMultiPions::NNInelasticIso(ener, 0)+CrossSectionsMultiPions::NNInelasticIso(ener, 2));
				if (xsinelas <= 1.e-9) return 0.;
			 G4double ratio=(NNToNNEta(particle1, particle2)-NNToNNEtaExclu(particle1, particle2))/xsinelas;
				if(s<6.25E6)
					return 0.;
				const G4double sigma = NNToNNEta(particle1, particle2) - NNToNNEtaExclu(particle1, particle2) - ratio*(NNToNNEtaOnePiOrDelta(particle1, particle2) + NNToNNEtaTwoPi(particle1, particle2) + NNToNNEtaThreePi(particle1, particle2));
				return ((sigma>1.e-9) ? sigma : 0.);
			}
			
			G4double CrossSectionsMultiPionsAndResonances::NNToNNEtaxPi(const G4int xpi, Particle const * const particle1, Particle const * const particle2) {
				//
				//     Nucleon-Nucleon producing one eta and xpi pions
				//
// assert(xpi>0 && xpi<=nMaxPiNN);
// assert(particle1->isNucleon() && particle2->isNucleon());
				
				const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 581.437; // 581.437 MeV translation to open pion production in NNEta
				if (ener < 2018.563) return 0.;
				const G4int i = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
    G4double xsinelas;
				if (i!=0) 
     xsinelas=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
    else 
					xsinelas=0.5*(CrossSectionsMultiPions::NNInelasticIso(ener, 0)+CrossSectionsMultiPions::NNInelasticIso(ener, 2));
				if (xsinelas <= 1.e-9) return 0.;
			 G4double ratio=(NNToNNEta(particle1, particle2)-NNToNNEtaExclu(particle1, particle2))/xsinelas;

				if (xpi == 1)
					return NNToNNEtaOnePi(particle1, particle2)*ratio;
				else if (xpi == 2)
					return NNToNNEtaTwoPi(particle1, particle2)*ratio;
				else if (xpi == 3)
					return NNToNNEtaThreePi(particle1, particle2)*ratio;
				else if (xpi == 4)
					return NNToNNEtaFourPi(particle1, particle2);
				else // should never reach this point
					return 0.;
			}

			
			G4double CrossSectionsMultiPionsAndResonances::NNToNDeltaEta(Particle const * const p1, Particle const * const p2) {
// assert(p1->isNucleon() && p2->isNucleon());
    const G4int i = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
				const G4double ener=KinematicsUtils::totalEnergyInCM(p1, p2) - 581.437; // 581.437 MeV translation to open pion production in NNEta
				if (ener < 2018.563) return 0.;
    G4double xsinelas;
				if (i!=0) 
     xsinelas=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
    else 
					xsinelas=0.5*(CrossSectionsMultiPions::NNInelasticIso(ener, 0)+CrossSectionsMultiPions::NNInelasticIso(ener, 2));
				if (xsinelas <= 1.e-9) return 0.;
			 G4double ratio=(NNToNNEta(p1, p2)-NNToNNEtaExclu(p1, p2))/xsinelas;
    G4double sigma = NNToNNEtaOnePiOrDelta(p1, p2)*ratio;
    if(i==0)
					sigma *= 0.5;
    return sigma;
			}

  
  G4double CrossSectionsMultiPionsAndResonances::NNToNNOmegaOnePi(Particle const * const particle1, Particle const * const particle2) {
   // Cross section for nucleon-nucleon producing one omega and one pion
   
   const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
   if (iso!=0)
    return 0.;
   
   const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 783.437; // 783.437 MeV translation to open pion production in NNOmega (= 2802. - 2018.563; 4074595.287720512986=2018.563*2018.563)
   if (ener < 2018.563) return 0.;
   
   const G4double xsiso2=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
   const G4double xsiso0=CrossSectionsMultiPions::NNInelasticIso(ener, 0);
   
   return 0.25*(CrossSectionsMultiPions::NNOnePiOrDelta(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNOnePiOrDelta(ener, 2, xsiso2));
  }
  
  G4double CrossSectionsMultiPionsAndResonances::NNToNNOmegaOnePiOrDelta(Particle const * const particle1, Particle const * const particle2) {
   const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 783.437; // 783.437 MeV translation to open pion production in NNOmega
   if (ener < 2018.563) return 0.;
   const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
   
   const G4double xsiso2=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
   if (iso != 0)
    return CrossSectionsMultiPions::NNOnePiOrDelta(ener, iso, xsiso2);
   else {
    const G4double xsiso0=CrossSectionsMultiPions::NNInelasticIso(ener, 0);
    return 0.5*(CrossSectionsMultiPions::NNOnePiOrDelta(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNOnePiOrDelta(ener, 2, xsiso2));
   }
  }
  
  G4double CrossSectionsMultiPionsAndResonances::NNToNNOmegaTwoPi(Particle const * const particle1, Particle const * const particle2) {
   //
   //     Nucleon-Nucleon producing one omega and two pions 
   //
   const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 783.437; // 783.437 MeV translation to open pion production in NNOmega
   if (ener < 2018.563) return 0.;
   const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
   
   
   const G4double xsiso2=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
   if (iso != 0) {
    return CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2);
   }
   else {
    const G4double xsiso0=CrossSectionsMultiPions::NNInelasticIso(ener, 0);
    return 0.5*(CrossSectionsMultiPions::NNTwoPi(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2));
   }
  }
  
  G4double CrossSectionsMultiPionsAndResonances::NNToNNOmegaThreePi(Particle const * const particle1, Particle const * const particle2) {
   //
   //     Nucleon-Nucleon producing one omega and three pions
   //
   
   const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 783.437; // 783.437 MeV translation to open pion production in NNOmega
   if (ener < 2018.563) return 0.;
   const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
   
   
   const G4double xsiso2=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
   const G4double xs1pi2=CrossSectionsMultiPions::NNOnePiOrDelta(ener, 2, xsiso2);
   const G4double xs2pi2=CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2);
   if (iso != 0)
    return CrossSectionsMultiPions::NNThreePi(ener, 2, xsiso2, xs1pi2, xs2pi2);
   else {
    const G4double xsiso0=CrossSectionsMultiPions::NNInelasticIso(ener, 0);
    const G4double xs1pi0=CrossSectionsMultiPions::NNOnePiOrDelta(ener, 0, xsiso0);
    const G4double xs2pi0=CrossSectionsMultiPions::NNTwoPi(ener, 0, xsiso0);
    return 0.5*(CrossSectionsMultiPions::NNThreePi(ener, 0, xsiso0, xs1pi0, xs2pi0)+ CrossSectionsMultiPions::NNThreePi(ener, 2, xsiso2, xs1pi2, xs2pi2));
   }
  }
  
  G4double CrossSectionsMultiPionsAndResonances::NNToNNOmegaFourPi(Particle const * const particle1, Particle const * const particle2) {
   //
   //     Nucleon-Nucleon producing one omega and four pions
   //
//jcd to be removed
//   return 0.;
//jcd    
   
   const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 783.437; // 783.437 MeV translation to open pion production in NNOmega
   if (ener < 2018.563) return 0.;
   const G4double s = ener*ener;
   const G4int i = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
   G4double xsinelas;
   if (i!=0) 
    xsinelas=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
   else 
    xsinelas=0.5*(CrossSectionsMultiPions::NNInelasticIso(ener, 0)+CrossSectionsMultiPions::NNInelasticIso(ener, 2));
   if (xsinelas <= 1.e-9) return 0.;
   G4double ratio=(NNToNNOmega(particle1, particle2)-NNToNNOmegaExclu(particle1, particle2))/xsinelas;
   if(s<6.25E6)
    return 0.;
   const G4double sigma = NNToNNOmega(particle1, particle2) - NNToNNOmegaExclu(particle1, particle2) - ratio*(NNToNNOmegaOnePiOrDelta(particle1, particle2) + NNToNNOmegaTwoPi(particle1, particle2) + NNToNNOmegaThreePi(particle1, particle2));
   return ((sigma>1.e-9) ? sigma : 0.);
  }
  
  G4double CrossSectionsMultiPionsAndResonances::NNToNNOmegaxPi(const G4int xpi, Particle const * const particle1, Particle const * const particle2) {
   //
   //     Nucleon-Nucleon producing one omega and xpi pions
   //
// assert(xpi>0 && xpi<=nMaxPiNN);
// assert(particle1->isNucleon() && particle2->isNucleon());
//jcd to be removed
//   return 0.;
//jcd    
   
   const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2) - 783.437; // 783.437 MeV translation to open pion production in NNOmega
   if (ener < 2018.563) return 0.;
   const G4int i = ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
   G4double xsinelas;
   if (i!=0) 
    xsinelas=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
   else 
    xsinelas=0.5*(CrossSectionsMultiPions::NNInelasticIso(ener, 0)+CrossSectionsMultiPions::NNInelasticIso(ener, 2));
   if (xsinelas <= 1.e-9) return 0.;
   G4double ratio=(NNToNNOmega(particle1, particle2)-NNToNNOmegaExclu(particle1, particle2))/xsinelas;
   
   if (xpi == 1)
    return NNToNNOmegaOnePi(particle1, particle2)*ratio;
   else if (xpi == 2)
    return NNToNNOmegaTwoPi(particle1, particle2)*ratio;
   else if (xpi == 3)
    return NNToNNOmegaThreePi(particle1, particle2)*ratio;
   else if (xpi == 4)
    return NNToNNOmegaFourPi(particle1, particle2);
   else // should never reach this point
    return 0.;
  }
  
  
  G4double CrossSectionsMultiPionsAndResonances::NNToNDeltaOmega(Particle const * const p1, Particle const * const p2) {
// assert(p1->isNucleon() && p2->isNucleon());
//jcd to be removed
//   return 0.;
//jcd    
   const G4int i = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
   const G4double ener=KinematicsUtils::totalEnergyInCM(p1, p2) - 783.437; // 783.437 MeV translation to open pion production in NNOmega
   if (ener < 2018.563) return 0.;
   G4double xsinelas;
   if (i!=0) 
    xsinelas=CrossSectionsMultiPions::NNInelasticIso(ener, 2);
   else 
    xsinelas=0.5*(CrossSectionsMultiPions::NNInelasticIso(ener, 0)+CrossSectionsMultiPions::NNInelasticIso(ener, 2));
   if (xsinelas <= 1.e-9) return 0.;
   G4double ratio=(NNToNNOmega(p1, p2)-NNToNNOmegaExclu(p1, p2))/xsinelas;
   G4double sigma = NNToNNOmegaOnePiOrDelta(p1, p2)*ratio;
   if(i==0)
    sigma *= 0.5;
   return sigma;
  }
  
	
} // namespace G4INCL

