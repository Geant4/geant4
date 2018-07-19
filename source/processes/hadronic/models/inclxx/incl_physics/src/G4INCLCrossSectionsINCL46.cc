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

#include "G4INCLCrossSectionsINCL46.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLLogger.hh"
// #include <cassert>

namespace G4INCL {

/*    G4double elasticNNHighEnergy(const G4double momentum) {
      return 77.0/(momentum + 1.5);
    }

    G4double elasticProtonNeutron(const G4double momentum) {
      if(momentum < 0.450) {
        const G4double alp = std::log(momentum);
        return 6.3555*std::exp(-3.2481*alp-0.377*alp*alp);
      } else if(momentum >= 0.450 && momentum < 0.8) {
        return (33.0 + 196.0 * std::sqrt(std::pow(std::abs(momentum - 0.95), 5)));
      } else if(momentum > 2.0) {
        return elasticNNHighEnergy(momentum);
      } else {
        return 31.0/std::sqrt(momentum);
      }
    }

    G4double elasticProtonProtonOrNeutronNeutron(const G4double momentum)
    {
      if(momentum < 0.440) {
        return 34.0*std::pow(momentum/0.4, -2.104);
      } else if(momentum < 0.8 && momentum >= 0.440) {
        return (23.5 + 1000.0*std::pow(momentum-0.7, 4));
      } else if(momentum < 2.0) {
        return (1250.0/(50.0 + momentum) - 4.0*std::pow(momentum-1.3, 2));
      } else {
        return elasticNNHighEnergy(momentum);
      }
    }

    G4double elasticNN(Particle const * const p1, Particle const * const p2) {
      G4double momentum = 0.0;
      momentum = 0.001 * KinematicsUtils::momentumInLab(p1, p2);
      if((p1->getType() == Proton && p2->getType() == Proton) ||
         (p1->getType() == Neutron && p2->getType() == Neutron)) {
        return elasticProtonProtonOrNeutronNeutron(momentum);
      } else if((p1->getType() == Proton && p2->getType() == Neutron) ||
                (p1->getType() == Neutron && p2->getType() == Proton)) {
        return elasticProtonNeutron(momentum);
      } else {
        INCL_ERROR("CrossSectionsINCL46::elasticNN: Bad input!" << '\n'
              << p1->print() << p2->print() << '\n');
      }
      return 0.0;
    }:*/

    G4double CrossSectionsINCL46::elasticNNLegacy(Particle const * const part1, Particle const * const part2) {


        G4int i = ParticleTable::getIsospin(part1->getType())
        + ParticleTable::getIsospin(part2->getType());

        /* The NN cross section is parametrised as a function of the lab momentum
         * of one of the nucleons. For NDelta or DeltaDelta, the physical
         * assumption is that the cross section is the same as NN *for the same
         * total CM energy*. Thus, we calculate s from the particles involved, and
         * we convert this value to the lab momentum of a nucleon *as if this were
         * an NN collision*.
         */

        const G4double s = KinematicsUtils::squareTotalEnergyInCM(part1, part2);
        G4double plab = 0.001*KinematicsUtils::momentumInLab(s, ParticleTable::effectiveNucleonMass, ParticleTable::effectiveNucleonMass);
        if(plab > 2.) { // NN, Delta-Nucleon and Delta-Delta for plab > 2.0 GeV
            return 77./(plab+1.5);
        }
        else if (part1->isNucleon() && part2->isNucleon()){ // NN
            if (i == 0) {   // pn
                if (plab < 0.450) {
                    G4double alp=std::log(plab);
                    return 6.3555*std::exp(-3.2481*alp-0.377*alp*alp);
                }
                else if (plab < 0.800) {
                    return (33.+196.*std::sqrt(std::pow(std::abs(plab-0.95),5)));
                }
                else {
                    return 31./std::sqrt(plab);
                }
            }
            else {   // nn and pp
                if (plab < 0.440) {
                    return 34.*std::pow(plab/0.4, (-2.104));
                }
                else if (plab < 0.800) {
                    return (23.5+1000.*std::pow(plab-0.7, 4));
                }
                else {
                    return (1250./(50.+plab)-4.*std::pow(plab-1.3, 2));
                }
            }
        }
        else {   // Delta-Nucleon and Delta-Delta
            if (plab < 0.440) {
                return 34.*std::pow(plab/0.4, (-2.104));
            }
            else if (plab < 0.800) {
                return (23.5+1000.*std::pow(plab-0.7, 4));
            }
            else {
                return (1250./(50.+plab)-4.*std::pow(plab-1.3, 2));
            }
        }
    }

  G4double CrossSectionsINCL46::deltaProduction(const G4int isospin, const G4double pLab) {
    G4double xs = 0.0;
// assert(isospin==-2 || isospin==0 || isospin==2);

    const G4double momentumGeV = 0.001 * pLab;
    if(pLab < 800.0) {
      return 0.0;
    }

    if(isospin==2 || isospin==-2) { // pp, nn
      if(pLab >= 2000.0) {
        xs = (41.0 + (60.0*momentumGeV - 54.0)*std::exp(-1.2*momentumGeV) - 77.0/(momentumGeV + 1.5));
      } else if(pLab >= 1500.0 && pLab < 2000.0) {
        xs = (41.0 + 60.0*(momentumGeV - 0.9)*std::exp(-1.2*momentumGeV) - 1250.0/(momentumGeV+50.0)+ 4.0*std::pow(momentumGeV - 1.3, 2));
      } else if(pLab < 1500.0) {
        xs = (23.5 + 24.6/(1.0 + std::exp(-10.0*momentumGeV + 12.0))
              -1250.0/(momentumGeV +50.0)+4.0*std::pow(momentumGeV - 1.3,2));
      }
    } else if(isospin==0) { // pn
      if(pLab >= 2000.0) {
        xs = (42.0 - 77.0/(momentumGeV + 1.5));
      } else if(pLab >= 1000.0 && pLab < 2000.0) {
        xs = (24.2 + 8.9*momentumGeV - 31.1/std::sqrt(momentumGeV));
      } else if(pLab < 1000.0) {
        xs = (33.0 + 196.0*std::sqrt(std::pow(std::abs(momentumGeV - 0.95),5))
              -31.1/std::sqrt(momentumGeV));
      }
    }

    if(xs < 0.0) return 0.0;
    else return xs;
  }

  G4double CrossSectionsINCL46::spnPiPlusPHE(const G4double x) {
    // HE and LE pi- p and pi+ n
    if(x <= 1750.0) {
      return -2.33730e-06*std::pow(x, 3)+1.13819e-02*std::pow(x,2)
        -1.83993e+01*x+9893.4;
    } else if(x > 1750.0 && x <= 2175.0) {
      return 1.13531e-06*std::pow(x, 3)-6.91694e-03*std::pow(x, 2)
        +1.39907e+01*x-9360.76;
    } else {
      return -3.18087*std::log(x)+52.9784;
    }
  }

  G4double CrossSectionsINCL46::spnPiMinusPHE(const G4double x) {
    // HE pi- p and pi+ n
    if(x <= 1475.0) {
      return 0.00120683*(x-1372.52)*(x-1372.52)+26.2058;
    } else if(x > 1475.0  && x <= 1565.0) {
      return 1.15873e-05*x*x+49965.6/((x-1519.59)*(x-1519.59)+2372.55);
    } else if(x > 1565.0 && x <= 2400.0) {
      return 34.0248+43262.2/((x-1681.65)*(x-1681.65)+1689.35);
    } else if(x > 2400.0 && x <= 7500.0) {
      return 3.3e-7*(x-7500.0)*(x-7500.0)+24.5;
    } else {
      return 24.5;
    }
  }

  G4double CrossSectionsINCL46::total(Particle const * const p1, Particle const * const p2) {
    G4double inelastic = 0.0;
    if(p1->isNucleon() && p2->isNucleon()) {
      inelastic = NNToNDelta(p1, p2);
    } else if((p1->isNucleon() && p2->isDelta()) ||
              (p1->isDelta() && p2->isNucleon())) {
      inelastic = NDeltaToNN(p1, p2);
    } else if((p1->isNucleon() && p2->isPion()) ||
              (p1->isPion() && p2->isNucleon())) {
      inelastic = piNToDelta(p1, p2);
    } else {
      inelastic = 0.0;
    }

    return inelastic + elastic(p1, p2);
  }

  G4double CrossSectionsINCL46::piNToDelta(Particle const * const particle1, Particle const * const particle2) {
    //      FUNCTION SPN(X,IND2T3,IPIT3,f17)
    // SIGMA(PI+ + P) IN THE (3,3) REGION
    // NEW FIT BY J.VANDERMEULEN  + FIT BY Th AOUST ABOVE (3,3) RES
    //                              CONST AT LOW AND VERY HIGH ENERGY
    //      COMMON/BL8/RATHR,RAMASS                                           REL21800
    //      integer f17
    // RATHR and RAMASS are always 0.0!!!

    G4double x = KinematicsUtils::totalEnergyInCM(particle1, particle2);
    if(x>10000.) return 0.0; // no cross section above this value

    G4int ipit3 = 0;
    G4int ind2t3 = 0;
    G4double ramass = 0.0;

    if(particle1->isPion()) {
      ipit3 = ParticleTable::getIsospin(particle1->getType());
    } else if(particle2->isPion()) {
      ipit3 = ParticleTable::getIsospin(particle2->getType());
    }

    if(particle1->isNucleon()) {
      ind2t3 = ParticleTable::getIsospin(particle1->getType());
    } else if(particle2->isNucleon()) {
      ind2t3 = ParticleTable::getIsospin(particle2->getType());
    }

    G4double y=x*x;
    G4double q2=(y-1076.0*1076.0)*(y-800.0*800.0)/y/4.0;
    if (q2 <= 0.) {
      return 0.0;
    }
    G4double q3 = std::pow(std::sqrt(q2),3);
    G4double f3 = q3/(q3 + 5832000.); // 5832000 = 180^3
    G4double spnResult = 326.5/(std::pow((x-1215.0-ramass)*2.0/(110.0-ramass), 2)+1.0);
    spnResult = spnResult*(1.0-5.0*ramass/1215.0);
    G4double cg = 4.0 + G4double(ind2t3)*G4double(ipit3);
    spnResult = spnResult*f3*cg/6.0;

    if(x < 1200.0  && spnResult < 5.0) {
      spnResult = 5.0;
    }

    // HE pi+ p and pi- n
    if(x > 1290.0) {
      if((ind2t3 == 1 && ipit3 == 2) || (ind2t3 == -1 && ipit3 == -2))
        spnResult=spnPiPlusPHE(x);
      else if((ind2t3 == 1 && ipit3 == -2) || (ind2t3 == -1 && ipit3 == 2))
        spnResult=spnPiMinusPHE(x);
      else if(ipit3 == 0) spnResult = (spnPiPlusPHE(x) + spnPiMinusPHE(x))/2.0; // (spnpipphe(x)+spnpimphe(x))/2.0
      else {
        INCL_ERROR("Unknown configuration!" << '\n');
      }
    }

    return spnResult;
  }

  G4double CrossSectionsINCL46::NDeltaToNN(Particle const * const p1, Particle const * const p2) {
    const G4int isospin = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
    if(isospin==4 || isospin==-4) return 0.0;

    G4double s = KinematicsUtils::squareTotalEnergyInCM(p1, p2);
    G4double Ecm = std::sqrt(s);
    G4int deltaIsospin;
    G4double deltaMass;
    if(p1->isDelta()) {
      deltaIsospin = ParticleTable::getIsospin(p1->getType());
      deltaMass = p1->getMass();
    } else {
      deltaIsospin = ParticleTable::getIsospin(p2->getType());
      deltaMass = p2->getMass();
    }

    if(Ecm <= 938.3 + deltaMass) {
      return 0.0;
    }

    if(Ecm < 938.3 + deltaMass + 2.0) {
      Ecm = 938.3 + deltaMass + 2.0;
      s = Ecm*Ecm;
    }

    const G4double x = (s - 4.*ParticleTable::effectiveNucleonMass2) /
      (s - std::pow(ParticleTable::effectiveNucleonMass + deltaMass, 2));
    const G4double y = s/(s - std::pow(deltaMass - ParticleTable::effectiveNucleonMass, 2));
    /* Concerning the way we calculate the lab momentum, see the considerations
     * in CrossSections::elasticNNLegacy().
     */
    const G4double pLab = KinematicsUtils::momentumInLab(s, ParticleTable::effectiveNucleonMass, ParticleTable::effectiveNucleonMass);
    G4double result = 0.5 * x * y * deltaProduction(isospin, pLab);
    result *= 3.*(32.0 + isospin * isospin * (deltaIsospin * deltaIsospin - 5))/64.0;
    result /= 1.0 + 0.25 * isospin * isospin;
    return result;
  }

  G4double CrossSectionsINCL46::NNToNDelta(Particle const * const p1, Particle const * const p2) {
// assert(p1->isNucleon() && p2->isNucleon());
    const G4double sqrts = KinematicsUtils::totalEnergyInCM(p1,p2);
    if(sqrts < ParticleTable::effectivePionMass + 2*ParticleTable::effectiveNucleonMass + 50.) { // approximately yields INCL4.6's hard-coded threshold in collis, 2065 MeV
      return 0.0;
    } else {
      const G4double pLab = KinematicsUtils::momentumInLab(p1,p2);
      const G4int isospin = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
      return deltaProduction(isospin, pLab);
    }
  }

  G4double CrossSectionsINCL46::elastic(Particle const * const p1, Particle const * const p2) {
//	  if(!p1->isPion() && !p2->isPion())
	  if((p1->isNucleon()||p1->isDelta()) && (p2->isNucleon()||p2->isDelta()))
      //    return elasticNN(p1, p2); // New implementation
      return elasticNNLegacy(p1, p2); // Translated from INCL4.6 FORTRAN
    else
      return 0.0; // No pion-nucleon elastic scattering
  }

  G4double CrossSectionsINCL46::calculateNNAngularSlope(G4double pl, G4int iso) {
    G4double x = 0.001 * pl; // Change to GeV
    if(iso != 0) {
      if(pl <= 2000.0) {
        x = std::pow(x, 8);
        return 5.5e-6 * x/(7.7 + x);
      } else {
        return (5.34 + 0.67*(x - 2.0)) * 1.0e-6;
      }
    } else {
      if(pl < 800.0) {
        G4double b = (7.16 - 1.63*x) * 1.0e-6;
        return b/(1.0 + std::exp(-(x - 0.45)/0.05));
      } else if(pl < 1100.0) {
        return (9.87 - 4.88 * x) * 1.0e-6;
      } else {
        return (3.68 + 0.76*x) * 1.0e-6;
      }
    }
    return 0.0; // Should never reach this point
  }


    G4double CrossSectionsINCL46::NNToxPiNN(const G4int, Particle const * const, Particle const * const) {
        return 0.;
    }
	
    G4double CrossSectionsINCL46::piNToxPiN(const G4int, Particle const * const, Particle const * const) {
        return 0.;
    }
	
    G4double CrossSectionsINCL46::piNToEtaN(Particle const * const, Particle const * const) {
		//
		//     Pion-Nucleon producing Eta cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsINCL46::piNToOmegaN(Particle const * const, Particle const * const) {
		//
		//     Pion-Nucleon producing Omega cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsINCL46::piNToEtaPrimeN(Particle const * const, Particle const * const) {
		//
		//     Pion-Nucleon producing EtaPrime cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsINCL46::etaNToPiN(Particle const * const, Particle const * const) {
		//
		//     Eta-Nucleon producing Pion cross sections
		//
		      return 0.;
    }

	
	   G4double CrossSectionsINCL46::etaNToPiPiN(Particle const * const, Particle const * const) {
		//
		//     Eta-Nucleon producing Two Pions cross sections
		//
		      return 0.;
	   }
	
    G4double CrossSectionsINCL46::omegaNToPiN(Particle const * const, Particle const * const) {
		//
		//     Omega-Nucleon producing Pion cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsINCL46::omegaNToPiPiN(Particle const * const, Particle const * const) {
		//
		//     Omega-Nucleon producing Two Pions cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsINCL46::etaPrimeNToPiN(Particle const * const, Particle const * const) {
		//
		//     EtaPrime-Nucleon producing Pion cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsINCL46::NNToNNEta(Particle const * const, Particle const * const) {
		//
		//     Nucleon-Nucleon producing Eta cross sections
		//
        return 0.;
    }
	
	   G4double CrossSectionsINCL46::NNToNNEtaExclu(Particle const * const, Particle const * const) {
		//
		//     Nucleon-Nucleon producing Eta cross sections
		//
							 return 0.;
	   }
	
   	G4double CrossSectionsINCL46::NNToNNEtaxPi(const G4int, Particle const * const, Particle const * const) {
	      	return 0.;
   	}

	   G4double CrossSectionsINCL46::NNToNDeltaEta(Particle const * const, Particle const * const) {
		//
		//     Nucleon-Nucleon producing N-Delta-Eta cross sections
		//
		      return 0.;
	}
 
    G4double CrossSectionsINCL46::NNToNNOmega(Particle const * const, Particle const * const) {
		//
		//     Nucleon-Nucleon producing Omega cross sections
		//
     return 0.;
    }
	
    G4double CrossSectionsINCL46::NNToNNOmegaExclu(Particle const * const, Particle const * const) {
		//
		//     Nucleon-Nucleon producing Omega cross sections
		//
     return 0.;
    }
	
    G4double CrossSectionsINCL46::NNToNNOmegaxPi(const G4int, Particle const * const, Particle const * const) {
     return 0.;
    }
 
    G4double CrossSectionsINCL46::NNToNDeltaOmega(Particle const * const, Particle const * const) {
  //
  //     Nucleon-Nucleon producing N-Delta-Omega cross sections
  //
     return 0.;
    }
	
	
	
    G4double CrossSectionsINCL46::NYelastic(Particle const * const , Particle const * const ) {
        //
        //      Hyperon-Nucleon elastic cross sections
        //
		return 0.;
    }

    G4double CrossSectionsINCL46::NKelastic(Particle const * const , Particle const * const ) {
        //
        //      Kaon-Nucleon elastic cross sections
        //
		return 0.;
	}

    G4double CrossSectionsINCL46::NKbelastic(Particle const * const , Particle const * const ) {
        //
        //      antiKaon-Nucleon elastic cross sections
        //
		return 0.;
	}
	
	G4double CrossSectionsINCL46::NNToNLK(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing N-Lambda-Kaon cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NNToNSK(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing N-Sigma-Kaon cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NNToNLKpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing N-Lambda-Kaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NNToNSKpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing N-Sigma-Kaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NNToNLK2pi(Particle const * const, Particle const * const) {
        //
        //     Nucleon-Nucleon producing N-Lambda-Kaon-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NNToNSK2pi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing N-Sigma-Kaon-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NNToNNKKb(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing Nucleon-Nucleon-Kaon-antiKaon cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NNToMissingStrangeness(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon missing strangeness production cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NDeltaToNLK(Particle const * const, Particle const * const) {
        // Nucleon-Delta producing Nucleon Lambda Kaon cross section
        return 0;
    }
    G4double CrossSectionsINCL46::NDeltaToNSK(Particle const * const, Particle const * const) {
        // Nucleon-Delta producing Nucleon Sigma Kaon cross section
        return 0;
    }
    G4double CrossSectionsINCL46::NDeltaToDeltaLK(Particle const * const, Particle const * const) {
        // Nucleon-Delta producing Delta Lambda Kaon cross section
        return 0;
    }
    G4double CrossSectionsINCL46::NDeltaToDeltaSK(Particle const * const, Particle const * const) {
        // Nucleon-Delta producing Delta Sigma Kaon cross section
        return 0;
    }
    
    G4double CrossSectionsINCL46::NDeltaToNNKKb(Particle const * const, Particle const * const) {
        // Nucleon-Delta producing Nucleon-Nucleon Kaon antiKaon cross section
        return 0;
    }

    G4double CrossSectionsINCL46::NpiToLK(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Lambda-Kaon cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NpiToSK(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Sigma-Kaon cross sections
        //
        return 0.;
    }
    G4double CrossSectionsINCL46::p_pimToSmKp(Particle const * const, Particle const * const) {
        return 0.;
    }
    G4double CrossSectionsINCL46::p_pimToSzKz(Particle const * const, Particle const * const) {
        return 0.;
    }
    G4double CrossSectionsINCL46::p_pizToSzKp(Particle const * const, Particle const * const) {
        return 0.;
    }

    G4double CrossSectionsINCL46::NpiToLKpi(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Lambda-Kaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NpiToSKpi(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Sigma-Kaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NpiToLK2pi(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Lambda-Kaon-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NpiToSK2pi(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Lambda-Kaon-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NpiToNKKb(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Nucleon-Kaon-antiKaon cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NpiToMissingStrangeness(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon missing strangeness production cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NLToNS(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Hyperon multiplet changing cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NSToNL(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Sigma quasi-elastic cross sections
        //
        return 0.;
    }
    
    G4double CrossSectionsINCL46::NSToNS(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Sigma quasi-elastic cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NKToNK(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Kaon quasi-elastic cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NKToNKpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Kaon producing Nucleon-Kaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NKToNK2pi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Kaon producing Nucleon-Kaon-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NKbToNKb(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon quasi-elastic cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NKbToSpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Sigma-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NKbToLpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Lambda-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NKbToS2pi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Sigma-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NKbToL2pi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Lambda-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NKbToNKbpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Nucleon-antiKaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsINCL46::NKbToNKb2pi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Nucleon-antiKaon-2pion cross sections
        //
        return 0.;
    }
	
} // namespace G4INCL

