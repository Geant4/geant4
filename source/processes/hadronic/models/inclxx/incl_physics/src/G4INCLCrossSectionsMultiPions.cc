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

#include "G4INCLCrossSectionsMultiPions.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLLogger.hh"
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

  const G4int CrossSectionsMultiPions::nMaxPiNN = 4;
  const G4int CrossSectionsMultiPions::nMaxPiPiN = 4;

  const G4double CrossSectionsMultiPions::s11pzOOT = 0.0035761542037692665889;
  const G4double CrossSectionsMultiPions::s01ppOOT = 0.003421025623481919853;
  const G4double CrossSectionsMultiPions::s01pzOOT = 0.0035739814152966403123;
  const G4double CrossSectionsMultiPions::s11pmOOT = 0.0034855350296270480281;
  const G4double CrossSectionsMultiPions::s12pmOOT = 0.0016672224074691565119;
  const G4double CrossSectionsMultiPions::s12ppOOT = 0.0016507643038726931312;
  const G4double CrossSectionsMultiPions::s12zzOOT = 0.0011111111111111111111;
  const G4double CrossSectionsMultiPions::s02pzOOT = 0.00125;
  const G4double CrossSectionsMultiPions::s02pmOOT = 0.0016661112962345883443;
  const G4double CrossSectionsMultiPions::s12mzOOT = 0.0017047391749062392793;

  CrossSectionsMultiPions::CrossSectionsMultiPions() :
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

  G4double CrossSectionsMultiPions::NNElastic(Particle const * const part1, Particle const * const part2) {

    /* The NN cross section is parametrised as a function of the lab momentum
     * of one of the nucleons. For NDelta or DeltaDelta, the physical
     * assumption is that the cross section is the same as NN *for the same
     * total CM energy*. Thus, we calculate s from the particles involved, and
     * we convert this value to the lab momentum of a nucleon *as if this were
     * an NN collision*.
     */
    const G4double s = KinematicsUtils::squareTotalEnergyInCM(part1, part2);

    if(part1->isNucleon() && part2->isNucleon()) {  // NN
      const G4int i = ParticleTable::getIsospin(part1->getType())
        + ParticleTable::getIsospin(part2->getType());
      return NNElasticFixed(s, i);
    }
    else {  // Nucleon-Delta and Delta-Delta
      const G4double plab = 0.001*KinematicsUtils::momentumInLab(s, ParticleTable::effectiveNucleonMass, ParticleTable::effectiveNucleonMass);
      if (plab < 0.440) {
        return 34.*std::pow(plab/0.4, (-2.104));
      }
      else if (plab < 0.800) {
        return 23.5+1000.*std::pow(plab-0.7, 4);
      }
      else if (plab <= 2.0) {
        return 1250./(50.+plab)-4.*std::pow(plab-1.3, 2);
      }
      else {
        return 77./(plab+1.5);
      }
    }
  }

    G4double CrossSectionsMultiPions::NNElasticFixed(const G4double s, const G4int i) {

      /* From NNElastic, with isospin fixed and for NN only.
      */

      G4double plab = 0.001*KinematicsUtils::momentumInLab(s, ParticleTable::effectiveNucleonMass, ParticleTable::effectiveNucleonMass);
      G4double sigma = 0.;

      if (i == 0) {  // pn
        if (plab < 0.446) {
          G4double alp=std::log(plab);
          sigma = 6.3555*std::exp(-3.2481*alp-0.377*alp*alp);
        }
        else if (plab < 0.851) {
          sigma = 33.+196.*std::pow(std::fabs(plab-0.95),2.5);
        }
        else if (plab <= 2.0) {
          sigma = 31./std::sqrt(plab);
        }
        else {
          sigma = 77./(plab+1.5);
        }
        //if(plab < 0.9 && plab > 0.802) sigma -= 0.1387*std::exp(-std::pow((plab-0.861),2)/0.0006861); //correction if totalcx-sumcx < 0.1
        //if(plab < 1.4 && plab > 1.31) sigma -= 0.1088*std::exp(-std::pow((plab-1.35),2)/0.00141); //correction if totalcx-sumcx < 0.1
        return sigma;
      }
      else {  // pp and nn
        if (plab < 0.440) {
          return 34.*std::pow(plab/0.4, (-2.104));
        }
        else if (plab < 0.8067) {
          return 23.5+1000.*std::pow(plab-0.7, 4);
        }
        else if (plab <= 2.0) {
          return 1250./(50.+plab)-4.*std::pow(plab-1.3, 2);
        }
        else if (plab <= 3.0956) {
          return 77./(plab+1.5);
        }
        else {
          G4double alp=std::log(plab);
          return 11.2+25.5*std::pow(plab, -1.12)+0.151*std::pow(alp, 2)-1.62*alp;
        }
      }
    }

    G4double CrossSectionsMultiPions::NNTot(Particle const * const part1, Particle const * const part2) {

        G4int i = ParticleTable::getIsospin(part1->getType())
        + ParticleTable::getIsospin(part2->getType());

        if(part1->isNucleon() && part2->isNucleon()) {  // NN
          const G4double s = KinematicsUtils::squareTotalEnergyInCM(part1, part2);
          return NNTotFixed(s, i);
        }
        else if (part1->isDelta() && part2->isDelta()) {  // Delta-Delta
            return elastic(part1, part2);
        }
        else {  // Nucleon-Delta
            return NDeltaToNN(part1, part2) + elastic(part1, part2);
        }
    }

    G4double CrossSectionsMultiPions::NNTotFixed(const G4double s, const G4int i) {

      /* From NNTot, with isospin fixed and for NN only.
      */

      G4double plab = 0.001*KinematicsUtils::momentumInLab(s, ParticleTable::effectiveNucleonMass, ParticleTable::effectiveNucleonMass);

      if (i == 0) {  // pn
        if (plab < 0.446) {
          G4double alp=std::log(plab);
          return 6.3555*std::exp(-3.2481*alp-0.377*std::pow(alp, 2));
        }
        else if (plab < 1.0) {
          return 33.+196.*std::sqrt(std::pow(std::fabs(plab-0.95),5));
        }
        else if (plab < 1.924) {
          return 24.2+8.9*plab;
        }
        else {
          G4double alp=std::log(plab);
          return 48.9-33.7*std::pow(plab, -3.08)+0.619*std::pow(alp, 2)-5.12*alp;
        }
      }
      else {  // pp and nn
        if (plab < 0.440) {
          return 34.*std::pow(plab/0.4, (-2.104));
        }
        else if (plab < 0.8734) {
          return 23.5+1000.*std::pow(plab-0.7, 4);
        }
        else if (plab < 1.5) {
          return 23.5+24.6/(1.+std::exp(-10.*(plab-1.2)));
        }
        else if (plab < 3.0044) {
          return 41.+60.*(plab-0.9)*std::exp(-1.2*plab);
        }
        else {
          G4double alp=std::log(plab);
          return 45.6+219.*std::pow(plab, -4.23)+0.41*std::pow(alp, 2)-3.41*alp;
        }
      }
    }

    G4double CrossSectionsMultiPions::NNInelasticIso(const G4double ener, const G4int iso) {

      const G4double s = ener*ener;
      G4double sincl;

      if (iso != 0) {
        if(s>=4074595.287720512986) { // plab>800 MeV/c
          sincl = NNTotFixed(s, 2)-NNElasticFixed(s, 2);
        }
        else {
          sincl =  0. ;
        }
      } else {
        if(s>=4074595.287720512986) { // plab>800 MeV/c
          sincl = 2*(NNTotFixed(s, 0)-NNElasticFixed(s, 0))-(NNTotFixed(s, 2)-NNElasticFixed(s, 2));
        }
        else {
          return 0. ;
        }
      }
      if (sincl < 0.) sincl = 0.;
      return sincl;
    }

    G4double CrossSectionsMultiPions::NNOnePiOrDelta(const G4double ener, const G4int iso, const G4double xsiso) {

        /* Article J. Physique 48 (1987)1901-1924 "Energy dependence of
         nucleon-cucleon inelastic total cross-sections."
         J. Bystricky, P. La France, F. Lehar, F. Perrot, T. Siemiarczuk & P. Winternitz
         S11PZ= section pp->pp pi0
         S01PP= section pp->pn pi+
         S01PZ= section pn->pn pi0
         S11PM= section pn->pp pi-
         S= X-Section, 1st number : 1 if pp and 0 if pn
         2nd number = number of pions, PP= pi+; PZ= pi0 ; PM= pi-
         */

        const G4double s = ener*ener;
        G4double plab = 0.001*KinematicsUtils::momentumInLab(s, ParticleTable::effectiveNucleonMass, ParticleTable::effectiveNucleonMass);

        G4double snnpit1=0.;
        G4double snnpit=0.;
        G4double s11pz=0.;
        G4double s01pp=0.;
        G4double s01pz=0.;
        G4double s11pm=0.;

        if ((iso != 0) && (plab < 2.1989)) {
            snnpit = xsiso - NNTwoPi(ener, iso, xsiso);
            if (snnpit < 1.e-8) snnpit=0.;
            return snnpit;
        }
        else if ((iso == 0) && (plab < 1.7369)) {
            snnpit = xsiso;
            if (snnpit < 1.e-8) snnpit=0.;
            return snnpit;
        }

//s11pz
        if (plab > 18.) {
            s11pz=55.185/std::pow((0.1412*plab+5),2);
        }
        else if (plab > 13.9) {
            G4double alp=std::log(plab);
            s11pz=6.67-13.3*std::pow(plab, -6.18)+0.456*alp*alp-3.29*alp;
        }
        else if (plab >= 0.7765) {
            const G4double b=BystrickyEvaluator<7>::eval(plab,s11pzOOT,s11pzHC);
            s11pz=b*b;
        }
//s01pp
        if (plab >= 0.79624) {
            const G4double b=BystrickyEvaluator<8>::eval(plab,s01ppOOT,s01ppHC);
            s01pp=b*b;
        }

// channel T=1
        snnpit1=s11pz+s01pp;
        if (snnpit1 < 1.e-8) snnpit1=0.;
        if (iso != 0) {
            return snnpit1;
        }

//s01pz
        if (plab > 4.5) {
            s01pz=15289.4/std::pow((11.573*plab+5),2);
        }
        else if (plab >= 0.777) {
            const G4double b=BystrickyEvaluator<4>::eval(plab,s01pzOOT,s01pzHC);
            s01pz=b*b;
        }
//s11pm
        if (plab > 14.) {
            s11pm=46.68/std::pow((0.2231*plab+5),2);
        }
        else if (plab >= 0.788) {
            const G4double b=BystrickyEvaluator<4>::eval(plab,s11pmOOT,s11pmHC);
            s11pm=b*b;
        }

// channel T=0
//        snnpit=s01pz+2*s11pm-snnpit1; //modif 2*(s01pz+2*s11pm)-snnpit1;
        snnpit = 2*(s01pz+2*s11pm)-snnpit1;
        if (snnpit < 1.e-8) snnpit=0.;
        return snnpit;
    }

    G4double CrossSectionsMultiPions::NNTwoPi(const G4double ener, const G4int iso, const G4double xsiso) {

        /* Article J. Physique 48 (1987)1901-1924 "Energy dependence of nucleon-cucleon inelastic total cross-sections."
           J. Bystricky, P. La France, F. Lehar, F. Perrot, T. Siemiarczuk & P. Winternitz
           S12PM : pp -> pp Pi+ Pi-
           S12ZZ : pp -> pp Pi0 Pi0
           S12PP : pp -> nn Pi+ Pi+
           S02PZ : pp -> pn Pi+ Pi0
           S02PM : pn -> pn Pi+ Pi-
           S12MZ : pn -> pp Pi- Pi0
        */

        const G4double s = ener*ener;
        G4double plab = 0.001*KinematicsUtils::momentumInLab(s, ParticleTable::effectiveNucleonMass, ParticleTable::effectiveNucleonMass);

        G4double snn2pit=0.;
        G4double s12pm=0.;
        G4double s12pp=0.;
        G4double s12zz=0.;
        G4double s02pz=0.;
        G4double s02pm=0.;
        G4double s12mz=0.;

        if (iso==0 && plab<3.33) {
            snn2pit = xsiso - NNOnePiOrDelta(ener, iso, xsiso);
            if (snn2pit < 1.e-8) snn2pit=0.;
            return snn2pit;
        }

        if (iso != 0) {
//s12pm
         if (plab > 15.) {
            s12pm=25.977/plab;
         }
         else if (plab >= 1.3817) {
            const G4double b=BystrickyEvaluator<5>::eval(plab,s12pmOOT,s12pmHC);
            s12pm=b*b;
         }
//s12pp
         if (plab > 10.) {
            s12pp=141.505/std::pow((-0.1016*plab-7),2);
         }
         else if (plab >= 1.5739) {
            const G4double b=BystrickyEvaluator<3>::eval(plab,s12ppOOT,s12ppHC);
            s12pp=b*b;
         }
        }
//s12zz
        if (plab > 4.) {
            s12zz=97.355/std::pow((1.1579*plab+5),2);
        }
        else if (plab >= 1.72207) {
            const G4double b=BystrickyEvaluator<4>::eval(plab,s12zzOOT,s12zzHC);
            s12zz=b*b;
        }
//s02pz
        if (plab > 4.5) {
            s02pz=178.082/std::pow((0.2014*plab+5),2);
        }
        else if (plab >= 1.5656) {
            const G4double b=BystrickyEvaluator<4>::eval(plab,s02pzOOT,s02pzHC);
            s02pz=b*b;
        }

// channel T=1
        if (iso != 0) {
            snn2pit=s12pm+s12pp+s12zz+s02pz;
            if (snn2pit < 1.e-8) snn2pit=0.;
            return snn2pit;
        }

//s02pm
        if (plab > 5.) {
            s02pm=135.826/std::pow(plab,2);
        }
        else if (plab >= 1.21925) {
            const G4double b=BystrickyEvaluator<6>::eval(plab,s02pmOOT,s02pmHC);
            s02pm=b*b;
        }
//s12mz
        if (plab >= 1.29269) {
            const G4double b=BystrickyEvaluator<4>::eval(plab,s12mzOOT,s12mzHC);
            s12mz=b*b;
        }

// channel T=0
//        snn2pit=3*(0.5*s02pm+0.5*s12mz-0.5*s02pz-s12zz);	//modif snn2pit=3*(s02pm+0.5*s12mz-0.5*s02pz-s12zz);
        snn2pit=3*(s02pm+0.5*s12mz-0.5*s02pz-s12zz);
        if (snn2pit < 1.e-8) snn2pit=0.;
        return snn2pit;
    }

    G4double CrossSectionsMultiPions::NNThreePi(const G4double ener, const G4int iso, const G4double xsiso, const G4double xs1pi, const G4double xs2pi) {

        const G4double s = ener*ener;
        G4double plab = 0.001*KinematicsUtils::momentumInLab(s, ParticleTable::effectiveNucleonMass, ParticleTable::effectiveNucleonMass);

        G4double snn3pit=0.;

        if (iso == 0) {
// channel T=0
            if (plab > 7.2355) {
                return 46.72/std::pow((plab - 5.8821),2);
            }
            else {
                snn3pit=xsiso-xs1pi-xs2pi;
                if (snn3pit < 1.e-8) snn3pit=0.;
                return snn3pit;
            }
        }
        else {
// channel T=1
            if (plab > 7.206) {
                return 5592.92/std::pow((plab+14.9764),2);
            }
            else if (plab > 2.1989){
                snn3pit=xsiso-xs1pi-xs2pi;
                if (snn3pit < 1.e-8) snn3pit=0.;
                return snn3pit;
            }
            else return snn3pit;
        }
    }

    G4double CrossSectionsMultiPions::NNOnePi(Particle const * const particle1, Particle const * const particle2) {
        // Cross section for nucleon-nucleon directly producing one pion

        const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());
        if (iso!=0) // If pp or nn we choose to always pass by the N-N to N-Delta channel
          return 0.;

        const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2);

        const G4double xsiso2=NNInelasticIso(ener, 2);
        const G4double xsiso0=NNInelasticIso(ener, 0);
        return 0.25*(NNOnePiOrDelta(ener, 0, xsiso0)+ NNOnePiOrDelta(ener, 2, xsiso2));
    }

    G4double CrossSectionsMultiPions::NNOnePiOrDelta(Particle const * const particle1, Particle const * const particle2) {
        // Cross section for nucleon-nucleon directly producing one pion or producing a nucleon-delta pair
        const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2);
        const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());

        const G4double xsiso2=NNInelasticIso(ener, 2);
        if (iso != 0)
          return NNOnePiOrDelta(ener, iso, xsiso2);
        else {
          const G4double xsiso0=NNInelasticIso(ener, 0);
          return 0.5*(NNOnePiOrDelta(ener, 0, xsiso0)+ NNOnePiOrDelta(ener, 2, xsiso2));
        }
    }

    G4double CrossSectionsMultiPions::NNTwoPi(Particle const * const particle1, Particle const * const particle2) {
        //
        //     Nucleon-Nucleon producing one pion cross sections
        //
        const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2);
        const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());


        const G4double xsiso2=NNInelasticIso(ener, 2);
        if (iso != 0) {
            return NNTwoPi(ener, 2, xsiso2);
        }
        else {
            const G4double xsiso0=NNInelasticIso(ener, 0);
            return 0.5*(NNTwoPi(ener, 0, xsiso0)+ NNTwoPi(ener, 2, xsiso2));
        }
        return 0.0; // Should never reach this point
    }

    G4double CrossSectionsMultiPions::NNThreePi(Particle const * const particle1, Particle const * const particle2) {
        //
        //     Nucleon-Nucleon producing one pion cross sections
        //

        const G4double ener=KinematicsUtils::totalEnergyInCM(particle1, particle2);
        const G4int iso=ParticleTable::getIsospin(particle1->getType()) + ParticleTable::getIsospin(particle2->getType());


        const G4double xsiso2=NNInelasticIso(ener, 2);
        const G4double xs1pi2=NNOnePiOrDelta(ener, 2, xsiso2);
        const G4double xs2pi2=NNTwoPi(ener, 2, xsiso2);
        if (iso != 0)
          return NNThreePi(ener, 2, xsiso2, xs1pi2, xs2pi2);
        else {
          const G4double xsiso0=NNInelasticIso(ener, 0);
          const G4double xs1pi0=NNOnePiOrDelta(ener, 0, xsiso0);
          const G4double xs2pi0=NNTwoPi(ener, 0, xsiso0);
          return 0.5*(NNThreePi(ener, 0, xsiso0, xs1pi0, xs2pi0)+ NNThreePi(ener, 2, xsiso2, xs1pi2, xs2pi2));
        }
    }

    G4double CrossSectionsMultiPions::NNFourPi(Particle const * const particle1, Particle const * const particle2) {
      const G4double s = KinematicsUtils::squareTotalEnergyInCM(particle1, particle2);
      if(s<6.25E6)
        return 0.;
      const G4double sigma = NNTot(particle1, particle2) - NNElastic(particle1, particle2) - NNOnePiOrDelta(particle1, particle2) - NNTwoPi(particle1, particle2) - NNThreePi(particle1, particle2);
      return ((sigma>1.e-9) ? sigma : 0.);
    }

    G4double CrossSectionsMultiPions::NNToxPiNN(const G4int xpi, Particle const * const particle1, Particle const * const particle2) {
      //
      //     Nucleon-Nucleon producing xpi pions cross sections
      //
// assert(xpi>0 && xpi<=nMaxPiNN);
// assert(particle1->isNucleon() && particle2->isNucleon());

      if (xpi == 1)
        return NNOnePi(particle1, particle2);
      else if (xpi == 2)
        return NNTwoPi(particle1, particle2);
      else if (xpi == 3)
        return NNThreePi(particle1, particle2);
      else if (xpi == 4)
        return NNFourPi(particle1, particle2);
      else // should never reach this point
        return 0.;
    }


  G4double CrossSectionsMultiPions::spnPiPlusPHE(const G4double x) {
    // HE and LE pi- p and pi+ n
    G4double ramass = 0.0;

    if(x <= 1306.78) {
       G4double y = x*x;
       G4double q2;
       q2=(y-std::pow(1076.0, 2))*(y-std::pow(800.0, 2))/(4.0*y);
       if (q2 > 0.) {
          G4double q3=std::pow(q2, 3./2.);
          G4double f3=q3/(q3+std::pow(180.0, 3));
	  G4double sdel;
	  sdel=326.5/(std::pow((x-1215.0-ramass)*2.0/110.0,2)+1.0);
	  return sdel*f3*(1.0-5.0*ramass/1215.0);
       }
       else {
          return 0;
       }
    }
    if(x <= 1754.0) {
      return -2.33730e-06*std::pow(x, 3)+1.13819e-02*std::pow(x,2)
        -1.83993e+01*x+9893.4;
    } else if (x <= 2150.0) {
      return 1.13531e-06*std::pow(x, 3)-6.91694e-03*std::pow(x, 2)
        +1.39907e+01*x-9360.76;
    } else {
      return -3.18087*std::log(x)+52.9784;
    }
  }

  G4double CrossSectionsMultiPions::spnPiMinusPHE(const G4double x) {
    // HE pi- p and pi+ n
    G4double ramass = 0.0;

    if(x <= 1275.8) {
       G4double y = x*x;
       G4double q2;
       q2=(y-std::pow(1076.0, 2))*(y-std::pow(800.0, 2))/(4.0*y);
       if (q2 > 0.) {
          G4double q3=std::pow(q2, 3./2.);
          G4double f3=q3/(q3+std::pow(180.0, 3));
	  G4double sdel;
	  sdel=326.5/(std::pow((x-1215.0-ramass)*2.0/110.0,2)+1.0);
	  return sdel*f3*(1.0-5.0*ramass/1215.0)/3.;
       }
       else {
          return 0;
       }
    }
    if(x <= 1495.0) {
      return 0.00120683*(x-1372.52)*(x-1372.52)+26.2058;
    } else if(x <= 1578.0) {
      return 1.15873e-05*x*x+49965.6/((x-1519.59)*(x-1519.59)+2372.55);
    } else if(x <= 2028.4) {
      return 34.0248+43262.2/((x-1681.65)*(x-1681.65)+1689.35);
    } else if(x <= 7500.0) {
      return 3.3e-7*(x-7500.0)*(x-7500.0)+24.5;
    } else {
      return 24.5;
    }
  }

  G4double CrossSectionsMultiPions::total(Particle const * const p1, Particle const * const p2) {
    G4double inelastic;
    if(p1->isNucleon() && p2->isNucleon()) {
      return NNTot(p1, p2);
    } else if((p1->isNucleon() && p2->isDelta()) ||
              (p1->isDelta() && p2->isNucleon())) {
      inelastic = NDeltaToNN(p1, p2);
    } else if((p1->isNucleon() && p2->isPion()) ||
              (p1->isPion() && p2->isNucleon())) {
      return piNTot(p1,p2);
    } else {
      inelastic = 0.;
    }

    return inelastic + elastic(p1, p2);
  }


	G4double CrossSectionsMultiPions::piNIne(Particle const * const particle1, Particle const * const particle2) {
		//      piN inelastic cross section (Delta excluded)
		
		const Particle *pion;
		const Particle *nucleon;
		if(particle1->isNucleon()) {
			nucleon = particle1;
			pion = particle2;
		} else {
			pion = particle1;
			nucleon = particle2;
		}
// assert(pion->isPion());
		
		const G4double pLab = KinematicsUtils::momentumInLab(pion, nucleon);
		
		// these limits correspond to sqrt(s)=1230 and 20000 MeV
		if(pLab>212677. || pLab<296.367)
			return 0.0;
		
		const G4int ipit3 = ParticleTable::getIsospin(pion->getType());
		const G4int ind2t3 = ParticleTable::getIsospin(nucleon->getType());
		const G4int cg = 4 + ind2t3*ipit3;
// assert(cg==2 || cg==4 || cg==6);
		
//		const G4double p1=1e-3*pLab;
//		const G4double p2=std::log(p1);
		G4double xpipp = 0.0;
		G4double xpimp = 0.0;
		
		if(cg!=2) {
			// x-section pi+ p inelastique :
			xpipp=piPluspIne(pion,nucleon);
			
			if(cg==6) // cas pi+ p et pi- n
				return xpipp;
		}
		
		// x-section pi- p inelastique :
		xpimp=piMinuspIne(pion,nucleon);
		
		if(cg==2) // cas pi- p et pi+ n
			return xpimp;
		else      // cas pi0 p et pi0 n
			return 0.5*(xpipp+xpimp);
	}
	
  G4double CrossSectionsMultiPions::piNToDelta(Particle const * const particle1, Particle const * const particle2) {
    //      piN Delta production

    G4double x = KinematicsUtils::totalEnergyInCM(particle1, particle2);
    if(x>20000.) return 0.0; // no cross section above this value

    G4int ipit3 = 0;
    G4int ind2t3 = 0;
    const G4double ramass = 0.0;

    if(particle1->isPion()) {
      ipit3 = ParticleTable::getIsospin(particle1->getType());
      ind2t3 = ParticleTable::getIsospin(particle2->getType());
    } else if(particle2->isPion()) {
      ipit3 = ParticleTable::getIsospin(particle2->getType());
      ind2t3 = ParticleTable::getIsospin(particle1->getType());
    }

    const G4double y=x*x;
    const G4double q2=(y-1076.0*1076.0)*(y-800.0*800.0)/y/4.0;
    if (q2 <= 0.) {
      return 0.0;
    }
    const G4double q3 = std::pow(std::sqrt(q2),3);
    const G4double f3 = q3/(q3 + 5832000.); // 5832000 = 180^3
    G4double sdelResult = 326.5/(std::pow((x-1215.0-ramass)*2.0/(110.0-ramass), 2)+1.0);
    sdelResult = sdelResult*(1.0-5.0*ramass/1215.0);
    const G4int cg = 4 + ind2t3*ipit3;
    sdelResult = sdelResult*f3*cg/6.0;

    return sdelResult;
  }

  G4double CrossSectionsMultiPions::piNTot(Particle const * const particle1, Particle const * const particle2) {
    //      FUNCTION SPN(X,IND2T3,IPIT3,f17)
    // SIGMA(PI+ + P) IN THE (3,3) REGION
    // NEW FIT BY J.VANDERMEULEN  + FIT BY Th AOUST ABOVE (3,3) RES
    //                              CONST AT LOW AND VERY HIGH ENERGY
    //      COMMON/BL8/RATHR,RAMASS                                           REL21800
    //      integer f17
    // RATHR and RAMASS are always 0.0!!!

    G4double x = KinematicsUtils::totalEnergyInCM(particle1, particle2);

    G4int ipit3 = 0;
    G4int ind2t3 = 0;

    if(particle1->isPion()) {
      ipit3 = ParticleTable::getIsospin(particle1->getType());
      ind2t3 = ParticleTable::getIsospin(particle2->getType());
    } else if(particle2->isPion()) {
      ipit3 = ParticleTable::getIsospin(particle2->getType());
      ind2t3 = ParticleTable::getIsospin(particle1->getType());
    }

    G4double spnResult=0.0;

    // HE pi+ p and pi- n
      if((ind2t3 == 1 && ipit3 == 2) || (ind2t3 == -1 && ipit3 == -2))
        spnResult=spnPiPlusPHE(x);
      else if((ind2t3 == 1 && ipit3 == -2) || (ind2t3 == -1 && ipit3 == 2))
        spnResult=spnPiMinusPHE(x);
      else if(ipit3 == 0) spnResult = (spnPiPlusPHE(x) + spnPiMinusPHE(x))/2.0; // (spnpipphe(x)+spnpimphe(x))/2.0
      else {
        INCL_ERROR("Unknown configuration!\n" << particle1->print() << particle2->print() << '\n');
      }

    return spnResult;
  }

  G4double CrossSectionsMultiPions::NDeltaToNN(Particle const * const p1, Particle const * const p2) {
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
    G4double sDelta;
    const G4double xsiso2=NNInelasticIso(Ecm, 2);
    if (isospin != 0)
      sDelta = NNOnePiOrDelta(Ecm, isospin, xsiso2);
    else {
      const G4double xsiso0=NNInelasticIso(Ecm, 0);
      sDelta = 0.25*(NNOnePiOrDelta(Ecm, 0, xsiso0)+ NNOnePiOrDelta(Ecm, 2, xsiso2));
    }
    G4double result = 0.5 * x * y * sDelta;
    /* modification for pion-induced cascade (see JC and MC LEMAIRE,NPA489(88)781
     * result=3.*result
     * pi absorption increased also for internal pions (7/3/01)
     */
    result *= 3.*(32.0 + isospin * isospin * (deltaIsospin * deltaIsospin - 5))/64.0;
    result /= 1.0 + 0.25 * (isospin * isospin);
    return result;
  }

  G4double CrossSectionsMultiPions::NNToNDelta(Particle const * const p1, Particle const * const p2) {
// assert(p1->isNucleon() && p2->isNucleon());
    const G4int isospin = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
    G4double sigma = NNOnePiOrDelta(p1, p2);
    if(isospin==0)
      sigma *= 0.5;
    return sigma;
  }

  G4double CrossSectionsMultiPions::elastic(Particle const * const p1, Particle const * const p2) {
//    if(!p1->isPion() && !p2->isPion()){
	  if((p1->isNucleon()||p1->isDelta()) && (p2->isNucleon()||p2->isDelta())){
      return NNElastic(p1, p2);
      }
//    else if (p1->isNucleon() || p2->isNucleon()){
	else if ((p1->isNucleon() && p2->isPion()) || (p2->isNucleon() && p1->isPion())){
      G4double pielas = piNTot(p1,p2) - piNIne(p1,p2) - piNToDelta(p1,p2);
        if (pielas < 0.){
            pielas = 0.;
        }
//        return piNTot(p1,p2) - piNIne(p1,p2) - piNToDelta(p1,p2);
        return pielas;
      }
    else {
       return 0.0;
      }
  }

  G4double CrossSectionsMultiPions::calculateNNAngularSlope(G4double pl, G4int iso) {
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


    G4double CrossSectionsMultiPions::piNToxPiN(const G4int xpi, Particle const * const particle1, Particle const * const particle2) {
        //
        //     pion-Nucleon producing xpi pions cross sections
        //
		const Particle *pion;
		const Particle *nucleon;
		if(particle1->isNucleon()) {
			nucleon = particle1;
			pion = particle2;
		} else {
			pion = particle1;
			nucleon = particle2;
		}
// assert(xpi>1 && xpi<=nMaxPiPiN);
// assert((particle1->isNucleon() && particle2->isPion()) || (particle1->isPion() && particle2->isNucleon()));
        const G4double plab = KinematicsUtils::momentumInLab(pion,nucleon);
		if (xpi == 2) {
			G4double OnePi=piNOnePi(particle1,particle2);
			if (OnePi < 1.e-09) OnePi = 0.;
            return OnePi;
        }
        else if (xpi == 3){
			G4double TwoPi=piNTwoPi(particle1,particle2);
			if (TwoPi < 1.e-09) TwoPi = 0.;									
            return TwoPi;
        }
        else if (xpi == 4) {
            G4double piNThreePi = piNIne(particle1,particle2) - piNOnePi(particle1,particle2) - piNTwoPi(particle1,particle2);
            if (piNThreePi < 1.e-09 || plab < 2000.) piNThreePi = 0.;									
            return piNThreePi;
        } else // should never reach this point
          return 0.0;
    }

	G4double CrossSectionsMultiPions::piNOnePi(Particle const * const particle1, Particle const * const particle2) {
		const Particle *pion;
		const Particle *nucleon;
		if(particle1->isNucleon()) {
			nucleon = particle1;
			pion = particle2;
		} else {
			pion = particle1;
			nucleon = particle2;
		}
// assert(pion->isPion());
		
		const G4double pLab = KinematicsUtils::momentumInLab(pion, nucleon);
		
		// this limit corresponds to sqrt(s)=1230 MeV
		if(pLab<296.367)
			return 0.0;
		
		const G4int ipi = ParticleTable::getIsospin(pion->getType());
		const G4int ind2 = ParticleTable::getIsospin(nucleon->getType());
		const G4int cg = 4 + ind2*ipi;
// assert(cg==2 || cg==4 || cg==6);
		
		//	const G4double p1=1e-3*pLab;
		G4double tamp6=0.;
		G4double tamp2=0.;
		const G4double elas = elastic(particle1, particle2);
		
		//   X-SECTION PI+ P INELASTIQUE :
		if(cg != 2) {
			tamp6=piPluspOnePi(particle1,particle2);
			if (cg == 6){ //   CAS PI+ P ET PI- N
				if(tamp6 >= elas && pLab < 410.) tamp6 = elas;
				return tamp6;
			}
		}
		
		//   X-SECTION PI- P INELASTIQUE :
		tamp2=piMinuspOnePi(particle1,particle2);
		if (tamp2 < 0.0) tamp2=0;
		
		if (cg == 2) //   CAS PI- P ET PI+ N
			return tamp2;
		else {       //   CAS PI0 P ET PI0 N
			G4double s1pin = 0.5*(tamp6+tamp2);
			const G4double inelastic = piNIne(particle1, particle2);
			if(s1pin >= elas && pLab < 410.) s1pin = 0.;
			if (s1pin > inelastic)
				s1pin = inelastic;
			return s1pin;
		}
	}
	
	G4double CrossSectionsMultiPions::piNTwoPi(Particle const * const particle1, Particle const * const particle2) {
		//
		//     pion-nucleon interaction, producing 2 pions
		//     fit from Landolt-Bornstein multiplied by factor determined with evaluation of total xs
		//
		
		const Particle *pion;
		const Particle *nucleon;
		if(particle1->isNucleon()) {
			nucleon = particle1;
			pion = particle2;
		} else {
			pion = particle1;
			nucleon = particle2;
		}
// assert(pion->isPion());
		
		const G4double pLab = KinematicsUtils::momentumInLab(pion, nucleon);
		const G4double elas = elastic(pion, nucleon);
		
		// this limit corresponds to sqrt(s)=1230 MeV
		if(pLab<296.367)
			return 0.0;
		
		const G4int ipi = ParticleTable::getIsospin(pion->getType());
		const G4int ind2 = ParticleTable::getIsospin(nucleon->getType());
		const G4int cg = 4 + ind2*ipi;
// assert(cg==2 || cg==4 || cg==6);
		
		G4double tamp6=0.;
		G4double tamp2=0.;
		
		//   X-SECTION PI+ P INELASTIQUE :
		if(cg!=2) {
			tamp6=piPluspTwoPi(particle1,particle2);
			if(cg==6){ //   CAS PI+ P ET PI- N
				if(tamp6 >= elas && pLab < 410.) tamp6 = 0.;
				return tamp6;}
		}
		
		//   X-SECTION PI- P INELASTIQUE :
		tamp2=piMinuspTwoPi(particle1,particle2);
		
		if(cg==2) //   CAS PI- P ET PI+ N
			return tamp2;
		else {    //   CAS PI0 P ET PI0 N
			const G4double s2pin=0.5*(tamp6+tamp2);
			return s2pin;
		}
	}
	
	G4double CrossSectionsMultiPions::piPluspIne(Particle const * const particle1, Particle const * const particle2) {
		//      piPlusP inelastic cross section (Delta excluded)
		
		const Particle *pion;
		const Particle *nucleon;
		if(particle1->isNucleon()) {
			nucleon = particle1;
			pion = particle2;
		} else {
			pion = particle1;
			nucleon = particle2;
		}
// assert(pion->isPion());
		
		const G4double pLab = KinematicsUtils::momentumInLab(pion, nucleon);
		
		// these limits correspond to sqrt(s)=1230 and 20000 MeV
		if(pLab>212677. || pLab<296.367)
			return 0.0;
		
//		const G4int ipit3 = ParticleTable::getIsospin(pion->getType());
//		const G4int ind2t3 = ParticleTable::getIsospin(nucleon->getType());
//		const G4int cg = 4 + ind2t3*ipit3;
//		assert(cg==2 || cg==4 || cg==6);
		
		const G4double p1=1e-3*pLab;
		const G4double p2=std::log(p1);
		G4double xpipp = 0.0;
		
		// x-section pi+ p inelastique :
		if(p1<=0.75)
			xpipp=17.965*std::pow(p1, 5.4606);
		else
			xpipp=24.3-12.3*std::pow(p1, -1.91)+0.324*p2*p2-2.44*p2;
		// cas pi+ p et pi- n
		return xpipp;
		
	}

	G4double CrossSectionsMultiPions::piMinuspIne(Particle const * const particle1, Particle const * const particle2) {
		//      piMinusp inelastic cross section (Delta excluded)
		
		const Particle *pion;
		const Particle *nucleon;
		if(particle1->isNucleon()) {
			nucleon = particle1;
			pion = particle2;
		} else {
			pion = particle1;
			nucleon = particle2;
		}
// assert(pion->isPion());
		
		const G4double pLab = KinematicsUtils::momentumInLab(pion, nucleon);
		
		// these limits correspond to sqrt(s)=1230 and 20000 MeV
		if(pLab>212677. || pLab<296.367)
			return 0.0;
		
//		const G4int ipit3 = ParticleTable::getIsospin(pion->getType());
//		const G4int ind2t3 = ParticleTable::getIsospin(nucleon->getType());
//		const G4int cg = 4 + ind2t3*ipit3;
//		assert(cg==2 || cg==4 || cg==6);
		
		const G4double p1=1e-3*pLab;
		const G4double p2=std::log(p1);
		G4double xpimp = 0.0;
		
		// x-section pi- p inelastique :
		if(p1 <= 0.4731)
			xpimp=0;
		else
			xpimp=26.6-7.18*std::pow(p1, -1.86)+0.327*p2*p2-2.81*p2;
		if(xpimp<0.)
			xpimp=0;
		
		// cas pi- p et pi+ n
		return xpimp;
		
	}

	G4double CrossSectionsMultiPions::piPluspOnePi(Particle const * const particle1, Particle const * const particle2) {
		const Particle *pion;
		const Particle *nucleon;
		if(particle1->isNucleon()) {
			nucleon = particle1;
			pion = particle2;
		} else {
			pion = particle1;
			nucleon = particle2;
		}
// assert(pion->isPion());
		
		const G4double pLab = KinematicsUtils::momentumInLab(pion, nucleon);
		
		// this limit corresponds to sqrt(s)=1230 MeV
		if(pLab<296.367)
			return 0.0;
		
		//	const G4int ipi = ParticleTable::getIsospin(pion->getType());
		//	const G4int ind2 = ParticleTable::getIsospin(nucleon->getType());
		//	const G4int cg = 4 + ind2*ipi;
		//	assert(cg==2 || cg==4 || cg==6);
		
		const G4double p1=1e-3*pLab;
		G4double tamp6=0.;
		
		//   X-SECTION PI+ P INELASTIQUE :
		if(pLab < 1532.52) // corresponds to sqrt(s)=1946 MeV
			tamp6=piPluspIne(particle1, particle2);
		else
			tamp6=0.204+18.2*std::pow(p1, -1.72)+6.33*std::pow(p1, -1.13);
		
		//   CAS PI+ P ET PI- N
		return tamp6;
		
	}

	G4double CrossSectionsMultiPions::piMinuspOnePi(Particle const * const particle1, Particle const * const particle2) {
		const Particle *pion;
		const Particle *nucleon;
		if(particle1->isNucleon()) {
			nucleon = particle1;
			pion = particle2;
		} else {
			pion = particle1;
			nucleon = particle2;
		}
// assert(pion->isPion());
		
		const G4double pLab = KinematicsUtils::momentumInLab(pion, nucleon);
		
		// this limit corresponds to sqrt(s)=1230 MeV
		if(pLab<296.367)
			return 0.0;
		
		//	const G4int ipi = ParticleTable::getIsospin(pion->getType());
		//	const G4int ind2 = ParticleTable::getIsospin(nucleon->getType());
		//	const G4int cg = 4 + ind2*ipi;
		//	assert(cg==2 || cg==4 || cg==6);
		
		const G4double p1=1e-3*pLab;
		G4double tamp2=0.;
		
		//   X-SECTION PI- P INELASTIQUE :
		if (pLab < 1228.06) // corresponds to sqrt(s)=1794 MeV
			tamp2=piMinuspIne(particle1, particle2);
		else
			tamp2=9.04*std::pow(p1, -1.17)+18.*std::pow(p1, -1.21); // tamp2=9.04*std::pow(p1, -1.17)+(13.5*std::pow(p1, -1.21))*4./3.;
		if (tamp2 < 0.0) tamp2=0;
		
		//   CAS PI- P ET PI+ N
		return tamp2;
	}

	G4double CrossSectionsMultiPions::piPluspTwoPi(Particle const * const particle1, Particle const * const particle2) {
		//
		//     pion-nucleon interaction, producing 2 pions
		//     fit from Landolt-Bornstein multiplied by factor determined with evaluation of total xs
		//
		
		const Particle *pion;
		const Particle *nucleon;
		if(particle1->isNucleon()) {
			nucleon = particle1;
			pion = particle2;
		} else {
			pion = particle1;
			nucleon = particle2;
		}
// assert(pion->isPion());
		
		const G4double pLab = KinematicsUtils::momentumInLab(pion, nucleon);
		
		// this limit corresponds to sqrt(s)=1230 MeV
		if(pLab<296.367)
			return 0.0;
		
		//	const G4int ipi = ParticleTable::getIsospin(pion->getType());
		//	const G4int ind2 = ParticleTable::getIsospin(nucleon->getType());
		//	const G4int cg = 4 + ind2*ipi;
		//	assert(cg==2 || cg==4 || cg==6);
		
		const G4double p1=1e-3*pLab;
		G4double tamp6=0.;
		
		//   X-SECTION PI+ P INELASTIQUE :
		if(pLab < 2444.7) // corresponds to sqrt(s)=2344 MeV
			tamp6=piPluspIne(particle1, particle2)-piPluspOnePi(particle1, particle2);
		else
			tamp6=1.59+25.5*std::pow(p1, -1.04); // tamp6=(0.636+10.2*std::pow(p1, -1.04))*15./6.;
		
		//   CAS PI+ P ET PI- N
		return tamp6;
	}
	
    G4double CrossSectionsMultiPions::piMinuspTwoPi(Particle const * const particle1, Particle const * const particle2) {
	//
	//     pion-nucleon interaction, producing 2 pions
	//     fit from Landolt-Bornstein multiplied by factor determined with evaluation of total xs
	//
	
	const Particle *pion;
	const Particle *nucleon;
	if(particle1->isNucleon()) {
		nucleon = particle1;
		pion = particle2;
	} else {
		pion = particle1;
		nucleon = particle2;
	}
// assert(pion->isPion());
	
	const G4double pLab = KinematicsUtils::momentumInLab(pion, nucleon);
	
	// this limit corresponds to sqrt(s)=1230 MeV
	if(pLab<296.367)
		return 0.0;
	
	//	const G4int ipi = ParticleTable::getIsospin(pion->getType());
	//	const G4int ind2 = ParticleTable::getIsospin(nucleon->getType());
	//	const G4int cg = 4 + ind2*ipi;
	//	assert(cg==2 || cg==4 || cg==6);
	
	const G4double p1=1e-3*pLab;
	G4double tamp2=0.;
	
	//   X-SECTION PI- P INELASTIQUE :
	if(pLab<2083.63) // corresponds to sqrt(s)=2195 MeV
		tamp2=piMinuspIne(particle1, particle2)-piMinuspOnePi(particle1, particle2);
	else
		tamp2=2.457794117647+18.066176470588*std::pow(p1, -0.92); // tamp2=(0.619+4.55*std::pow(p1, -0.92))*135./34.;
	
	//   CAS PI- P ET PI+ N
	return tamp2;
}



	
    G4double CrossSectionsMultiPions::piNToEtaN(Particle const * const, Particle const * const) {
		//
		//     Pion-Nucleon producing Eta cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsMultiPions::piNToOmegaN(Particle const * const, Particle const * const) {
		//
		//     Pion-Nucleon producing Omega cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsMultiPions::piNToEtaPrimeN(Particle const * const, Particle const * const) {
		//
		//     Pion-Nucleon producing EtaPrime cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsMultiPions::etaNToPiN(Particle const * const, Particle const * const) {
		//
		//     Eta-Nucleon producing Pion cross sections
		//
		      return 0.;
    }

	
	   G4double CrossSectionsMultiPions::etaNToPiPiN(Particle const * const, Particle const * const) {
		//
		//     Eta-Nucleon producing Two Pions cross sections
		//
		      return 0.;
	   }
	
	
    G4double CrossSectionsMultiPions::omegaNToPiN(Particle const * const, Particle const * const) {
		//
		//     Omega-Nucleon producing Pion cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsMultiPions::omegaNToPiPiN(Particle const * const, Particle const * const) {
		//
		//     Omega-Nucleon producing Two Pions cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsMultiPions::etaPrimeNToPiN(Particle const * const, Particle const * const) {
		//
		//     EtaPrime-Nucleon producing Pion cross sections
		//
        return 0.;
    }
	
    G4double CrossSectionsMultiPions::NNToNNEta(Particle const * const, Particle const * const) {
		//
		//     Nucleon-Nucleon producing Eta cross sections
		//
        return 0.;
    }
	
	G4double CrossSectionsMultiPions::NNToNNEtaExclu(Particle const * const, Particle const * const) {
		//
		//     Nucleon-Nucleon producing Eta cross sections
		//
	    return 0.;
	   }
	
	G4double CrossSectionsMultiPions::NNToNNEtaxPi(const G4int, Particle const * const, Particle const * const) {
	    return 0.;
	   }

   	G4double CrossSectionsMultiPions::NNToNDeltaEta(Particle const * const, Particle const * const) {
		//
		//     Nucleon-Nucleon producing N-Delta-Eta cross sections
		//
		return 0.;
		}

    G4double CrossSectionsMultiPions::NNToNNOmega(Particle const * const, Particle const * const) {
		//
		//     Nucleon-Nucleon producing Omega cross sections
		//
     return 0.;
    }
	
    G4double CrossSectionsMultiPions::NNToNNOmegaExclu(Particle const * const, Particle const * const) {
		//
		//     Nucleon-Nucleon producing Omega cross sections
		//
     return 0.;
    }
	
    G4double CrossSectionsMultiPions::NNToNNOmegaxPi(const G4int, Particle const * const, Particle const * const) {
     return 0.;
    }
 
    G4double CrossSectionsMultiPions::NNToNDeltaOmega(Particle const * const, Particle const * const) {
  //
  //     Nucleon-Nucleon producing N-Delta-Omega cross sections
  //
     return 0.;
    }




    G4double CrossSectionsMultiPions::NYelastic(Particle const * const , Particle const * const ) {
        //
        //      Hyperon-Nucleon elastic cross sections
        //
		return 0.;
    }

    G4double CrossSectionsMultiPions::NKelastic(Particle const * const , Particle const * const ) {
        //
        //      Kaon-Nucleon elastic cross sections
        //
		return 0.;
	}

    G4double CrossSectionsMultiPions::NKbelastic(Particle const * const , Particle const * const ) {
        //
        //      antiKaon-Nucleon elastic cross sections
        //
		return 0.;
	}


	G4double CrossSectionsMultiPions::NNToNLK(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing N-Lambda-Kaon cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NNToNSK(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing N-Sigma-Kaon cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NNToNLKpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing N-Lambda-Kaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NNToNSKpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing N-Sigma-Kaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NNToNLK2pi(Particle const * const, Particle const * const) {
        //
        //     Nucleon-Nucleon producing N-Lambda-Kaon-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NNToNSK2pi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing N-Sigma-Kaon-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NNToNNKKb(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon producing Nucleon-Nucleon-Kaon-antiKaon cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NNToMissingStrangeness(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Nucleon missing strangeness production cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NDeltaToNLK(Particle const * const, Particle const * const) {
        // Nucleon-Delta producing Nucleon Lambda Kaon cross section
        return 0;
    }
    G4double CrossSectionsMultiPions::NDeltaToNSK(Particle const * const, Particle const * const) {
        // Nucleon-Delta producing Nucleon Sigma Kaon cross section
        return 0;
    }
    G4double CrossSectionsMultiPions::NDeltaToDeltaLK(Particle const * const, Particle const * const) {
        // Nucleon-Delta producing Delta Lambda Kaon cross section
        return 0;
    }
    G4double CrossSectionsMultiPions::NDeltaToDeltaSK(Particle const * const, Particle const * const) {
        // Nucleon-Delta producing Delta Sigma Kaon cross section
        return 0;
    }
    
    G4double CrossSectionsMultiPions::NDeltaToNNKKb(Particle const * const, Particle const * const) {
        // Nucleon-Delta producing Nucleon-Nucleon Kaon antiKaon cross section
        return 0;
    }


    G4double CrossSectionsMultiPions::NpiToLK(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Lambda-Kaon cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NpiToSK(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Sigma-Kaon cross sections
        //
        return 0.;
    }
    G4double CrossSectionsMultiPions::p_pimToSmKp(Particle const * const, Particle const * const) {
        return 0.;
    }
    G4double CrossSectionsMultiPions::p_pimToSzKz(Particle const * const, Particle const * const) {
        return 0.;
    }
    G4double CrossSectionsMultiPions::p_pizToSzKp(Particle const * const, Particle const * const) {
        return 0.;
    }

    G4double CrossSectionsMultiPions::NpiToLKpi(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Lambda-Kaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NpiToSKpi(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Sigma-Kaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NpiToLK2pi(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Lambda-Kaon-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NpiToSK2pi(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Lambda-Kaon-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NpiToNKKb(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon producing Nucleon-Kaon-antiKaon cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NpiToMissingStrangeness(Particle const * const, Particle const * const) {
        //
        //      Pion-Nucleon missing strangeness production cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NLToNS(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Hyperon multiplet changing cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NSToNL(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Sigma quasi-elastic cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NSToNS(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Sigma quasi-elastic cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NKToNK(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Kaon quasi-elastic cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NKToNKpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Kaon producing Nucleon-Kaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NKToNK2pi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-Kaon producing Nucleon-Kaon-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NKbToNKb(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon quasi-elastic cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NKbToSpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Sigma-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NKbToLpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Lambda-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NKbToS2pi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Sigma-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NKbToL2pi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Lambda-2pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NKbToNKbpi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Nucleon-antiKaon-pion cross sections
        //
        return 0.;
    }

    G4double CrossSectionsMultiPions::NKbToNKb2pi(Particle const * const, Particle const * const) {
        //
        //      Nucleon-antiKaon producing Nucleon-antiKaon-2pion cross sections
        //
        return 0.;
    }



	
} // namespace G4INCL

