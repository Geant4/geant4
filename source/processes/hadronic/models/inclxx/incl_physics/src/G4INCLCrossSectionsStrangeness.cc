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

/** \file G4INCLCrossSectionsStrangeness.cc
 * \brief Multipion, mesonic Resonances and strange cross sections
 *
 * \date 1st March 2016
 * \author Jason Hirtz
 */

#include "G4INCLCrossSectionsStrangeness.hh"
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
    
    const G4int CrossSectionsStrangeness::nMaxPiNN = 4;
    const G4int CrossSectionsStrangeness::nMaxPiPiN = 4;

    CrossSectionsStrangeness::CrossSectionsStrangeness() :
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

	/// \brief redefining previous cross sections

    G4double CrossSectionsStrangeness::total(Particle const * const p1, Particle const * const p2) {
        G4double inelastic;
        if(p1->isNucleon() && p2->isNucleon()) {
            return CrossSectionsMultiPions::NNTot(p1, p2);
        } else if((p1->isNucleon() && p2->isDelta()) ||
                  (p1->isDelta() && p2->isNucleon())) {
            inelastic = CrossSectionsMultiPions::NDeltaToNN(p1, p2) + NDeltaToNLK(p1, p2) + NDeltaToNSK(p1, p2) + NDeltaToDeltaLK(p1, p2) + NDeltaToDeltaSK(p1, p2) + NDeltaToNNKKb(p1, p2);
        } else if((p1->isNucleon() && p2->isPion()) ||
                  (p1->isPion() && p2->isNucleon())) {
            return CrossSectionsMultiPions::piNTot(p1,p2);
        } else if((p1->isNucleon() && p2->isEta()) ||
                  (p1->isEta() && p2->isNucleon())) {
            inelastic = CrossSectionsMultiPionsAndResonances::etaNToPiN(p1,p2) + CrossSectionsMultiPionsAndResonances::etaNToPiPiN(p1,p2);
        } else if((p1->isNucleon() && p2->isOmega()) ||
                  (p1->isOmega() && p2->isNucleon())) {
            inelastic = CrossSectionsMultiPionsAndResonances::omegaNInelastic(p1,p2);
        } else if((p1->isNucleon() && p2->isEtaPrime()) ||
                  (p1->isEtaPrime() && p2->isNucleon())) {
            inelastic = CrossSectionsMultiPionsAndResonances::etaPrimeNToPiN(p1,p2);
        } else if((p1->isNucleon() && p2->isLambda()) ||
                  (p1->isLambda() && p2->isNucleon())) {
            inelastic = NLToNS(p1,p2);
        } else if((p1->isNucleon() && p2->isSigma()) ||
                  (p1->isSigma() && p2->isNucleon())) {
            inelastic = NSToNL(p1,p2) + NSToNS(p1,p2);
        } else if((p1->isNucleon() && p2->isKaon()) ||
                  (p1->isKaon() && p2->isNucleon())) {
            inelastic = NKToNK(p1,p2) + NKToNKpi(p1,p2) + NKToNK2pi(p1,p2);
        } else if((p1->isNucleon() && p2->isAntiKaon()) ||
                  (p1->isAntiKaon() && p2->isNucleon())) {
            inelastic = NKbToLpi(p1,p2) + NKbToSpi(p1,p2) + NKbToL2pi(p1,p2) + NKbToS2pi(p1,p2) + NKbToNKb(p1,p2) + NKbToNKbpi(p1,p2) + NKbToNKb2pi(p1,p2);
        } else {
            inelastic = 0.;
        }
        return inelastic + elastic(p1, p2);
    }    

    G4double CrossSectionsStrangeness::elastic(Particle const * const p1, Particle const * const p2) {
        if((p1->isNucleon()||p1->isDelta()) && (p2->isNucleon()||p2->isDelta())){ // N-N, N-Delta, Delta-Delta
            return CrossSectionsMultiPions::elastic(p1, p2);
        }
        else if ((p1->isNucleon() && p2->isPion()) || (p2->isNucleon() && p1->isPion())){
            return CrossSectionsMultiPions::elastic(p1, p2);
        }
        else if ((p1->isNucleon() && p2->isEta()) || (p2->isNucleon() && p1->isEta())){
            return CrossSectionsMultiPionsAndResonances::etaNElastic(p1, p2);
        }
        else if ((p1->isNucleon() && p2->isHyperon()) || (p2->isNucleon() && p1->isHyperon())){
            return NYelastic(p1, p2);
        }
        else if ((p1->isNucleon() && p2->isKaon()) || (p2->isNucleon() && p1->isKaon())){
            return NKelastic(p1, p2);
        }
        else if ((p1->isNucleon() && p2->isAntiKaon()) || (p2->isNucleon() && p1->isAntiKaon())){
            return NKbelastic(p1, p2);
        }
        else {
            return 0.0;
        }
    }

    G4double CrossSectionsStrangeness::piNToxPiN(const G4int xpi, Particle const * const particle1, Particle const * const particle2) {
        //
        //     pion-Nucleon producing xpi pions cross sections (corrected due to new particles)
        //
// assert(xpi>1 && xpi<=nMaxPiPiN);
// assert((particle1->isNucleon() && particle2->isPion()) || (particle1->isPion() && particle2->isNucleon()));
                
            const G4double oldXS2Pi=CrossSectionsMultiPions::piNToxPiN(2,particle1, particle2);
        const G4double oldXS3Pi=CrossSectionsMultiPions::piNToxPiN(3,particle1, particle2);
        const G4double oldXS4Pi=CrossSectionsMultiPions::piNToxPiN(4,particle1, particle2);
        const G4double xsEta=CrossSectionsMultiPionsAndResonances::piNToEtaN(particle1, particle2);
        const G4double xsOmega=CrossSectionsMultiPionsAndResonances::piNToOmegaN(particle1, particle2);
        const G4double xs1=NpiToLK(particle2, particle1);
        const G4double xs2=NpiToSK(particle1, particle2);
        const G4double xs3=NpiToLKpi(particle1, particle2);
        const G4double xs4=NpiToSKpi(particle1, particle2);
        const G4double xs5=NpiToLK2pi(particle1, particle2);
        const G4double xs6=NpiToSK2pi(particle1, particle2);
        const G4double xs7=NpiToNKKb(particle1, particle2);
        const G4double xs8=NpiToMissingStrangeness(particle1, particle2);
        const G4double xs0 = xs1 + xs2 + xs3 + xs4 + xs5 + xs6 + xs7 +xs8;
        G4double newXS2Pi=0.;    
        G4double newXS3Pi=0.;    
        G4double newXS4Pi=0.;
        
        if (xpi == 2) {
            if (oldXS4Pi != 0.)
                newXS2Pi=oldXS2Pi;
            else if (oldXS3Pi != 0.) {
                newXS3Pi=oldXS3Pi-xsEta-xsOmega-xs0;
                if (newXS3Pi < 1.e-09)
                    newXS2Pi=oldXS2Pi-(xsEta+xsOmega+xs0-oldXS3Pi);
                else
                    newXS2Pi=oldXS2Pi;
            }
            else { 
                newXS2Pi=oldXS2Pi-xsEta-xsOmega-xs0;
                if (newXS2Pi < 1.e-09 && newXS2Pi!=0.){
                    newXS2Pi=0.;
                }
            }
            return newXS2Pi;
        }
        else if (xpi == 3) {
            if (oldXS4Pi != 0.) {
                newXS4Pi=oldXS4Pi-xsEta-xsOmega-xs0;
                if (newXS4Pi < 1.e-09)
                    newXS3Pi=oldXS3Pi-(xsEta+xsOmega+xs0-oldXS4Pi);
                else
                    newXS3Pi=oldXS3Pi;
            }
            else { 
                newXS3Pi=oldXS3Pi-xsEta-xsOmega-xs0;
                if (newXS3Pi < 1.e-09){
                    newXS3Pi=0.;
                }
            }
            return newXS3Pi;
        }
        else if (xpi == 4) {
            newXS4Pi=oldXS4Pi-xsEta-xsOmega-xs0;
            if (newXS4Pi < 1.e-09){
                    newXS4Pi=0.;
                }
            return newXS4Pi;
        }
        else // should never reach this point
            return 0.;
    }

    G4double CrossSectionsStrangeness::NNToxPiNN(const G4int xpi, Particle const * const particle1, Particle const * const particle2) {
        //
        //     Nucleon-Nucleon producing xpi pions cross sections corrected
        //
// assert(xpi>0 && xpi<=nMaxPiNN);
// assert(particle1->isNucleon() && particle2->isNucleon());
        
        G4double oldXS1Pi=CrossSectionsMultiPions::NNToxPiNN(1,particle1, particle2);
        G4double oldXS2Pi=CrossSectionsMultiPions::NNToxPiNN(2,particle1, particle2);
        G4double oldXS3Pi=CrossSectionsMultiPions::NNToxPiNN(3,particle1, particle2);
        G4double oldXS4Pi=CrossSectionsMultiPions::NNToxPiNN(4,particle1, particle2);
        G4double xsEta=CrossSectionsMultiPionsAndResonances::NNToNNEta(particle1, particle2)+CrossSectionsMultiPionsAndResonances::NNToNNOmega(particle1, particle2);;
        
        const G4double xs1=NNToNLK(particle1, particle2);
        const G4double xs2=NNToNSK(particle1, particle2);
        const G4double xs3=NNToNLKpi(particle1, particle2);
        const G4double xs4=NNToNSKpi(particle1, particle2);
        const G4double xs5=NNToNLK2pi(particle1, particle2);
        const G4double xs6=NNToNSK2pi(particle1, particle2);
        const G4double xs7=NNToNNKKb(particle1, particle2);
        const G4double xs8=NNToMissingStrangeness(particle1, particle2);
        const G4double xs0 = xs1 + xs2 + xs3 + xs4 + xs5 + xs6 + xs7 + xs8;
        G4double newXS1Pi=0.;    
        G4double newXS2Pi=0.;    
        G4double newXS3Pi=0.;    
        G4double newXS4Pi=0.;    
            
        if (xpi == 1) {
            if (oldXS4Pi != 0. || oldXS3Pi != 0.)
                newXS1Pi=oldXS1Pi;
            else if (oldXS2Pi != 0.) { 
                newXS2Pi=oldXS2Pi-xsEta-xs0;
                if (newXS2Pi < 0.)
                    newXS1Pi=oldXS1Pi-(xsEta+xs0-oldXS2Pi);
                else
                    newXS1Pi=oldXS1Pi;
            }
            else 
                newXS1Pi=oldXS1Pi-xsEta-xs0;
            return newXS1Pi;
        }
        else if (xpi == 2) {
            if (oldXS4Pi != 0.){
                newXS2Pi=oldXS2Pi;
            }
            else if (oldXS3Pi != 0.) {
                newXS3Pi=oldXS3Pi-xsEta-xs0;
                if (newXS3Pi < 0.)
                    newXS2Pi = oldXS2Pi-(xsEta+xs0-oldXS3Pi);
                else
                    newXS2Pi = oldXS2Pi;
            }
            else { 
                newXS2Pi = oldXS2Pi-xsEta-xs0;
                if (newXS2Pi < 0.)
                    newXS2Pi = 0.;
            }
            return newXS2Pi;
        }
        else if (xpi == 3) {
            if (oldXS4Pi != 0.) {
                newXS4Pi=oldXS4Pi-xsEta-xs0;
             if (newXS4Pi < 0.)
                 newXS3Pi=oldXS3Pi-(xsEta+xs0-oldXS4Pi);
             else
                 newXS3Pi=oldXS3Pi;
            }
            else { 
                newXS3Pi=oldXS3Pi-xsEta-xs0;                
                if (newXS3Pi < 0.)
                    newXS3Pi=0.;
            }
            return newXS3Pi;
        }
        else if (xpi == 4) {
            newXS4Pi=oldXS4Pi-xsEta-xs0;
            if (newXS4Pi < 0.)
                newXS4Pi=0.;
        return newXS4Pi;
        }

        else // should never reach this point
            return 0.;
    }

	/// \brief elastic cross sections

    G4double CrossSectionsStrangeness::NYelastic(Particle const * const p1, Particle const * const p2) {
        //
        //      Hyperon-Nucleon elastic cross sections
        //
// assert((p1->isNucleon() && p2->isHyperon()) || (p1->isHyperon() && p2->isNucleon()));
        
        G4double sigma = 0.;        
        
        const Particle *hyperon;
        const Particle *nucleon;
        
        if (p1->isHyperon()) {
            hyperon = p1;
            nucleon = p2;
        }
        else {
            hyperon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = KinematicsUtils::momentumInLab(hyperon, nucleon); // MeV
        
        if (pLab < 145.)
            sigma = 200.;
        else if (pLab < 425.)
            sigma = 869.*std::exp(-pLab/100.);
        else if (pLab < 30000.)
            sigma = 12.8*std::exp(-6.2e-5*pLab);
        else    
            sigma=0.;
            
        if (sigma < 0.) sigma = 0.; // should never happen
        return sigma;
    }

    G4double CrossSectionsStrangeness::NKelastic(Particle const * const p1, Particle const * const p2) {
        //
        //      Kaon-Nucleon elastic cross sections
        //
// assert((p1->isNucleon() && p2->isKaon()) || (p1->isKaon() && p2->isNucleon()));
        
        G4double sigma=0.;        
        
        const Particle *kaon;
        const Particle *nucleon;
        
        if (p1->isKaon()) {
            kaon = p1;
            nucleon = p2;
        }
        else {
            kaon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = KinematicsUtils::momentumInLab(kaon, nucleon); // MeV
        
        if (pLab < 935.)
            sigma = 12.;
        else if (pLab < 2080.)
            sigma = 17.4-3.*std::exp(6.3e-4*pLab);
        else if (pLab < 5500.)
            sigma = 832.*std::pow(pLab,-0.64);
        else if (pLab < 30000.)
            sigma = 3.36;
        else    
            sigma=0.;
            
        if (sigma < 0.) sigma = 0.; // should never happen
        return sigma;
    }

    G4double CrossSectionsStrangeness::NKbelastic(Particle const * const p1, Particle const * const p2) {
        //
        //      antiKaon-Nucleon elastic cross sections
        //
// assert((p1->isNucleon() && p2->isAntiKaon()) || (p1->isAntiKaon() && p2->isNucleon()));
        
        G4double sigma=0.;        
        
        const Particle *antikaon;
        const Particle *nucleon;
        
        if (p1->isAntiKaon()) {
            antikaon = p1;
            nucleon = p2;
        }
        else {
            antikaon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antikaon, nucleon); // GeV
        
        if(pLab > 1E-6) // sigma = 287.823 [mb] -> rise very slowly, not cut needed
           sigma = 6.132*std::pow(pLab,-0.2437)+12.98*std::exp(-std::pow(pLab-0.9902,2)/0.05558)+2.928*std::exp(-std::pow(pLab-1.649,2)/0.772)+564.3*std::exp(-std::pow(pLab+0.9901,2)/0.5995);
            
        if (sigma < 0.) sigma = 0.;  // should never happen
        return sigma;
    }
    
	/// \brief NN to strange cross sections

    G4double CrossSectionsStrangeness::NNToNLK(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Nucleon producing N-Lambda-Kaon cross sections
        //
        // ratio
        // p p (1)    p n (1)
        //
        // p p -> p L K+ (1)
        // p n -> p L K0 (1/2)
        // p n -> n L K+ (1/2)
// assert(p1->isNucleon() && p2->isNucleon());
        
        const Particle *particle1;
        const Particle *particle2;
        
        if(p2->getType() == Proton && p1->getType() == Neutron){
            particle1 = p2;
            particle2 = p1;
        }
        else{
            particle1 = p1;
            particle2 = p2;
        }
        
        G4double sigma = 0.;
        
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(particle1, particle2); // GeV
        
        if(particle2->getType() == Proton){
            if(pLab < 2.3393) return 0.;
            else if (pLab < 30.) sigma = 1.11875*std::pow((pLab-2.3393),1.0951)/std::pow((pLab+2.3393),2.0958); // pLab = 30 Gev -> exess of energie = 5 GeV
            else return 0.;
        }
        else{
            if(pLab < 2.3508) return 0.;
            else if (pLab < 30.) sigma = 1.11875*std::pow((pLab-2.3508),1.0951)/std::pow((pLab+2.3508),2.0958);
            else return 0.;
        }
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NNToNSK(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Nucleon producing N-Sigma-Kaon cross sections
        //
        // Meson symmetry
        // pp->pS+K0 (1/4)
        // pp->pS0K+ (1/8) // HEM
        // pp->pS0K+ (1/4) // Data
        // pp->nS+K+ (1)
        //
        // pn->nS+K0 (1/4)
        // pn->pS-K+ (1/4)
        // pn->nS0K+ (5/8)
        // pn->pS0K0 (5/8)
        //
// assert(p1->isNucleon() && p2->isNucleon());
        
        const Particle *particle1;
        const Particle *particle2;
        
        if(p2->getType() == Proton && p1->getType() == Neutron){
            particle1 = p2;
            particle2 = p1;
        }
        else{
            particle1 = p1;
            particle2 = p2;
        }
        
        G4double sigma = 0.;
        
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(particle1, particle2); // GeV
        
        if(pLab < 2.593)
            return 0.;
        
        if(p2->getType() == p1->getType())
//            sigma = 1.375*2*6.38*std::pow(pLab-2.57,2.1)/std::pow(pLab,4.162);
            sigma = 1.5*6.38*std::pow(pLab-2.593,2.1)/std::pow(pLab,4.162);
        else
            sigma = 1.75*6.38*std::pow(pLab-2.593,2.1)/std::pow(pLab,4.162);
            
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NNToNLKpi(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Nucleon producing N-Lambda-Kaon-pion cross sections
        //
        // ratio (pure NN -> DLK)
        // pp (12)    pn (8)
        //
        // pp -> p pi+ L K0 (9)(3)
        // pp -> p pi0 L K+ (2)(1*2/3)
        // pp -> n pi+ L K+ (1)(1*1/3)
        //
        // pn -> p pi- L K+ (2)(2*1/3)
        // pn -> n pi0 L K+ (4)(2*2/3)
        // pn -> p pi0 L K0 (4)
        // pn -> n pi+ L K0 (2)
        
// assert(p1->isNucleon() && p2->isNucleon());
                
        G4double sigma = 0.;
        G4double ratio = 0.;
        G4double ratio1 = 0.;
        G4double ratio2 = 0.;
        const G4double ener = KinematicsUtils::totalEnergyInCM(p1, p2) - 540.;
        if( ener < p1->getMass() + p2->getMass())
            return 0;
        const G4int iso = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const G4double xsiso2 = CrossSectionsMultiPions::NNInelasticIso(ener, 2);
        if (iso != 0){
            ratio1 = CrossSectionsMultiPions::NNOnePiOrDelta(ener, iso, xsiso2);
            ratio2 = CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2);
        }
        else {
            const G4double xsiso0 = CrossSectionsMultiPions::NNInelasticIso(ener, 0);
            ratio1 = 0.5*(CrossSectionsMultiPions::NNOnePiOrDelta(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNOnePiOrDelta(ener, 2, xsiso2));
            ratio2 = 0.5*(CrossSectionsMultiPions::NNTwoPi(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2));
        }
        
        if( ratio1 == 0 || ratio2 == 0)
            return 0.;
        
        ratio = ratio2/ratio1;
        
        sigma = ratio * NNToNLK(p1,p2) * 3;
        
/*        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(p1, p2); // GeV
        if(pLab <= 2.77) return 0.;
        sigma = 0.4 * std::pow(pLab-2.77,1.603)/std::pow(pLab,1.492);*/
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NNToNSKpi(Particle const * const p1, Particle const * const p2) {
        //
        // Nucleon-Nucleon producing N-Sigma-Kaon-pion cross sections
        //
        // ratio (pure NN -> DSK)
        // pp (36)    pn (36)
        //
        // pp -> p pi+ S- K+ (9)
        // pp -> p pi+ S0 K0 (9)
        // pp -> p pi0 S+ K0 (4)
        // pp -> n pi+ S+ K0 (2)
        // pp -> p pi0 S0 K+ (4)
        // pp -> n pi+ S0 K+ (2)
        // pp -> p pi- S+ K+ (2)
        // pp -> n pi0 S+ K+ (4)
        
        // pn -> p pi0 S- K+ (4)
        // pn -> n pi+ S- K+ (2)
        // pn -> p pi0 S0 K0 (2)
        // pn -> n pi+ S0 K0 (1)
        // pn -> p pi+ S- K0 (9)
        
// assert(p1->isNucleon() && p2->isNucleon());
                
        G4double sigma = 0.;
        G4double ratio = 0.;
        G4double ratio1 = 0.;
        G4double ratio2 = 0.;
        const G4double ener = KinematicsUtils::totalEnergyInCM(p1, p2) - 620.;
        if( ener < p1->getMass() + p2->getMass())
            return 0;
        const G4int iso = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const G4double xsiso2 = CrossSectionsMultiPions::NNInelasticIso(ener, 2);
        if (iso != 0){
            ratio1 = CrossSectionsMultiPions::NNOnePiOrDelta(ener, iso, xsiso2);
            ratio2 = CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2);
        }
        else {
            const G4double xsiso0 = CrossSectionsMultiPions::NNInelasticIso(ener, 0);
            ratio1 = 0.5*(CrossSectionsMultiPions::NNOnePiOrDelta(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNOnePiOrDelta(ener, 2, xsiso2));
            ratio2 = 0.5*(CrossSectionsMultiPions::NNTwoPi(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2));
        }
        
        if( ratio1 == 0 || ratio2 == 0)
            return 0.;
        
        ratio = ratio2/ratio1;
        
        sigma = ratio * NNToNSK(p1,p2) * 3;
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NNToNLK2pi(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Nucleon producing N-Lambda-Kaon-pion cross sections
        //
// assert(p1->isNucleon() && p2->isNucleon());
                
        G4double sigma = 0.;
        G4double ratio = 0.;
        G4double ratio1 = 0.;
        G4double ratio2 = 0.;
        const G4double ener = KinematicsUtils::totalEnergyInCM(p1, p2) - 675.;
        if( ener < p1->getMass() + p2->getMass())
            return 0;
        const G4int iso = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const G4double xsiso2 = CrossSectionsMultiPions::NNInelasticIso(ener, 2);
        if (iso != 0){
            ratio1 = CrossSectionsMultiPions::NNOnePiOrDelta(ener, iso, xsiso2);
            ratio2 = CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2);
        }
        else {
            const G4double xsiso0 = CrossSectionsMultiPions::NNInelasticIso(ener, 0);
            ratio1 = 0.5*(CrossSectionsMultiPions::NNOnePiOrDelta(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNOnePiOrDelta(ener, 2, xsiso2));
            ratio2 = 0.5*(CrossSectionsMultiPions::NNTwoPi(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2));
        }
        
        if( ratio1 == 0 || ratio2 == 0)
            return 0.;
        
        ratio = ratio2/ratio1;
        
        sigma = ratio * NNToNLKpi(p1,p2);
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NNToNSK2pi(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Nucleon producing N-Sigma-Kaon-pion cross sections
        //
// assert(p1->isNucleon() && p2->isNucleon());
                
        G4double sigma = 0.;
        G4double ratio = 0.;
        G4double ratio1 = 0.;
        G4double ratio2 = 0.;
        const G4double ener = KinematicsUtils::totalEnergyInCM(p1, p2) - 755.;
        if( ener < p1->getMass() + p2->getMass())
            return 0;
        const G4int iso = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const G4double xsiso2 = CrossSectionsMultiPions::NNInelasticIso(ener, 2);
        if (iso != 0){
            ratio1 = CrossSectionsMultiPions::NNOnePiOrDelta(ener, iso, xsiso2);
            ratio2 = CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2);
        }
        else {
            const G4double xsiso0 = CrossSectionsMultiPions::NNInelasticIso(ener, 0);
            ratio1 = 0.5*(CrossSectionsMultiPions::NNOnePiOrDelta(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNOnePiOrDelta(ener, 2, xsiso2));
            ratio2 = 0.5*(CrossSectionsMultiPions::NNTwoPi(ener, 0, xsiso0)+ CrossSectionsMultiPions::NNTwoPi(ener, 2, xsiso2));
        }
        
        if( ratio1 == 0 || ratio2 == 0)
            return 0.;
        
        ratio = ratio2/ratio1;
        
        sigma = ratio * NNToNSKpi(p1,p2);
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NNToNNKKb(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Nucleon producing Nucleon-Nucleon-Kaon-antiKaon cross sections
        //
        // Channel strongly resonant; fit from Sibirtesev - Z. Phys. A 358, 101-106 (1997) (eq.21)
        // ratio
        // pp (6) pn (13)*2
        // pp -> pp K+ K- (1)
        // pp -> pp K0 K0 (1)
        // pp -> pn K+ K0 (4)
        // pn -> pp K0 K- (4)
        // pn -> pn K+ K- (9)
        //
        
// assert(p1->isNucleon() && p2->isNucleon());
                
        G4double sigma = 0.;
        const G4int iso = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        const G4double ener = 0.001*KinematicsUtils::totalEnergyInCM(p1, p2); // GeV
        
        if(ener < 2.872)
            return 0.;
        
        if(iso == 0)
            sigma = 26 * 5./19. * 0.3 *std::pow(1.-2.872*2.872/(ener*ener),3.)*std::pow(2.872*2.872/(ener*ener),0.8);
        else
            sigma = 6 * 5./19. * 0.3 *std::pow(1.-2.872*2.872/(ener*ener),3.)*std::pow(2.872*2.872/(ener*ener),0.8);
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NNToMissingStrangeness(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Nucleon missing strangeness production cross sections
        //
// assert(p1->isNucleon() && p2->isNucleon());
        
        G4double sigma = 0.;
        
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(p1,p2); // GeV
        const G4int iso = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        if(pLab < 6.) return 0.;
        
        if(iso == 0){
            if(pLab < 30.) sigma = 10.15*std::pow((pLab - 6.),2.157)/std::pow(pLab,2.333);
            else return 0.;
        }
        else{
            if(pLab < 30.) sigma = 8.12*std::pow((pLab - 6.),2.157)/std::pow(pLab,2.333);
            else return 0.;
        }
        return sigma;
    }

	/** \brief NDelta to strange cross sections
	 * 
	 * No experimental data
	 * Parametrization from Phys.Rev.C 59 1 (369) (1999)
	 * 
	 * Correction are applied on the isospin symetry provided in the paper
	 */
	 
    G4double CrossSectionsStrangeness::NDeltaToNLK(Particle const * const p1, Particle const * const p2) {
        // Nucleon-Delta producing Nucleon Lambda Kaon cross section
        //
        // XS from K. Tsushima, A. Sibirtsev, A. W. Thomas, and G. Q. Li. Phys.Rev.C 59, 369
        // 
        // D++ n -> p L K+ (3)
        //
        // D+  p -> p L K+ (1)
        //
        // D+  n -> p L K0 (1)
        // D+  n -> n L K+ (1)
        //return 0.;
// assert((p1->isNucleon() && p2->isResonance()) || (p2->isNucleon() && p1->isResonance()));
        
        const G4int iso = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        if(std::abs(iso) == 4) return 0.;
                
        G4double sigma = 0.;
        
        const G4double s = KinematicsUtils::squareTotalEnergyInCM(p1 ,p2); // MeV^2
        const G4double s0 = 6.511E6; // MeV^2
        
        if(s <= s0) return 0.;
        
        sigma = 4.*4.169*std::pow(s/s0-1,2.227)*std::pow(s0/s,2.511);
        
        //const G4double pLab = sdt::sqrt(s*s/(4*ParticleTable::effectiveNucleonMass2)-s)*0.001;
        //sigma = 3*1.11875*std::pow((pLab-2.3508),1.0951)/std::pow((pLab+2.3508),2.0958); // NDelta sim to NN
        
        if(iso == 0){ // D+  n
            sigma *= 2./6.;
        }
        else if (ParticleTable::getIsospin(p1->getType()) == ParticleTable::getIsospin(p2->getType())){// D+  p
            sigma *= 1./6.;
        }
        else{// D++ n
            sigma *= 3./6.;
        }
        return sigma;
    }
    G4double CrossSectionsStrangeness::NDeltaToNSK(Particle const * const p1, Particle const * const p2) {
        // Nucleon-Delta producing Nucleon Sigma Kaon cross section
        //
        // XS from K. Tsushima, A. Sibirtsev, A. W. Thomas, and G. Q. Li. Phys.Rev.C 59, 369 ( X 1.25 (124/99) for isospin consideration)
        //
        // D++ p -> p S+ K+ (6)
        //
        // D++ n -> p S+ K0 (3) ****
        // D++ n -> p S0 K+ (3)
        // D++ n -> n S+ K+ (3)
        //
        // D+  p -> p S+ K0 (2)
        // D+  p -> p S0 K+ (2)
        // D+  p -> n S+ K+ (3)
        //
        // D+  n -> p S0 K0 (3)
        // D+  n -> p S- K+ (2)
        // D+  n -> n S+ K0 (2)
        // D+  n -> n S0 K+ (2)
        
// assert((p1->isNucleon() && p2->isResonance()) || (p2->isNucleon() && p1->isResonance()));
        
        G4double sigma = 0.;
        
        const G4double s = KinematicsUtils::squareTotalEnergyInCM(p1 ,p2); // Mev^^2
        const G4double s0 = 6.935E6; // Mev^2
        const G4int iso = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        if(s <= s0)
            return 0.;
        
        sigma = 11.*39.54*std::pow(s/s0-1,2.799)*std::pow(s0/s,6.303);
        
        //const G4double pLab = sdt::sqrt(s*s/(4*ParticleTable::effectiveNucleonMass2)-s)*0.001;
        //sigma = 22./12./2. * 4.75*6.38*std::pow(pLab-2.593,2.1)/std::pow(pLab,4.162); // NDelta sim to NN
        
        if(iso == 0)// D+  n
            sigma *= 9./31.;
        else if (ParticleTable::getIsospin(p1->getType()) == ParticleTable::getIsospin(p2->getType()))// D+  p
            sigma *= 7./31.;
        else if (std::abs(iso) == 2)// D++ n
            sigma *= 9./31.;
        else // D++ p
            sigma *= 6./31.;
        
        return sigma;
    }
    G4double CrossSectionsStrangeness::NDeltaToDeltaLK(Particle const * const p1, Particle const * const p2) {
        // Nucleon-Delta producing Delta Lambda Kaon cross section
        //
        // XS from K. Tsushima, A. Sibirtsev, A. W. Thomas, and G. Q. Li. Phys.Rev.C 59, 369
        //
        // ratio
        // D++ p -> L K+ D++ (4)
        //
        // D++ n -> L K+ D+  (3)
        // D++ n -> L K0 D++ (4)
        //
        // D+  p -> L K0 D++ (3)
        // D+  p -> L K+ D+  (2)
        //
        // D+  n -> L K+ D0  (4)
        // D+  n -> L K0 D+  (2)
        //return 0.;
        
// assert((p1->isNucleon() && p2->isResonance()) || (p2->isNucleon() && p1->isResonance()));
        
        const G4double s = KinematicsUtils::squareTotalEnergyInCM(p1 ,p2); // Mev^^2
        const G4double s0 = 8.096E6; // Mev^2
        const G4int iso = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        if(s <= s0)
            return 0.;
        
        G4double sigma = 7.*2.679*std::pow(s/s0-1,2.280)*std::pow(s0/s,5.086);
        
        if(iso == 0)// D+  n
            sigma *= 6./22.;
        else if (ParticleTable::getIsospin(p1->getType()) == ParticleTable::getIsospin(p2->getType()))// D+  p
            sigma *=  5./22.;
        else if (std::abs(iso) == 2)// D++ n
            sigma *=  7./22.;
        else // D++ p
            sigma *=  4./22.;
        
        return sigma;
    }
    G4double CrossSectionsStrangeness::NDeltaToDeltaSK(Particle const * const p1, Particle const * const p2) {
        // Nucleon-Delta producing Delta Sigma Kaon cross section
        //
        // XS from K. Tsushima, A. Sibirtsev, A. W. Thomas, and G. Q. Li. Phys.Rev.C 59, 369
        //
        // D++ p (9)
        // D++ n (15)
        // D+  p (11)
        // D+  n (13)
        //
        // ratio
        // D++ p -> S+ K+ D+  (a)	(2)
        // D++ p -> S0 K+ D++ (b)	(1)
        // D++ p -> S+ K0 D++ (c)	(6)
        //
        // D++ n -> S+ K+ D0 *(d)*	(2)
        // D++ n -> S0 K+ D+  (e)	(4)
        // D++ n -> S- K+ D++ (f)	(6)(c)*
        // D++ n -> S+ K0 D+  (a)*	(2)
        // D++ n -> S0 K0 D++ (b)*	(1)*
        //
        // D+  p -> S+ K+ D0  (i)	(2)*
        // D+  p -> S0 K+ D+  (j)	(1)
        // D+  p -> S- K+ D++ (k)	(2)(g=a)*
        // D+  p -> S+ K0 D+  (l)	(2)
        // D+  p -> S0 K0 D++ (m)	(4)(e)*
        //
        // D+  n -> S+ K+ D- *(d)*	(2)
        // D+  n -> S0 K+ D0  (o)	(4)
        // D+  n -> S- K+ D+  (p)	(2)*
        // D+  n -> S+ K0 D0  (i)*	(2)*
        // D+  n -> S0 K0 D+  (j)*	(1)*
        // D+  n -> S- K0 D++ (k)*	(2)*
        //return 0.;
        
// assert((p1->isNucleon() && p2->isResonance()) || (p2->isNucleon() && p1->isResonance()));
        
        const G4double s = KinematicsUtils::squareTotalEnergyInCM(p1 ,p2); // Mev^^2
        const G4double s0 = 8.568E6; // Mev^2
        const G4int iso = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        if(s <= s0)
            return 0.;
        
        G4double sigma = 19.*21.18*std::pow(s/s0-1,2.743)*std::pow(s0/s,8.407);
        
        if(iso == 0)// D+  n
            sigma *= 13./48.;
        else if (ParticleTable::getIsospin(p1->getType()) == ParticleTable::getIsospin(p2->getType()))// D+  p
            sigma *=  11./48.;
        else if (std::abs(iso) == 2)// D++ n
            sigma *=  15./48.;
        else // D++ p
            sigma *=  9./48.;
        
        return sigma;
    }
    
    G4double CrossSectionsStrangeness::NDeltaToNNKKb(Particle const * const p1, Particle const * const p2) {
        // Nucleon-Delta producing Nucleon-Nucleon Kaon antiKaon cross section
        //
        // Total = sigma(NN->NNKKb)*10
        //
        // D++ p (6)
        // D++ n (9)
        // D+  p (7)
        // D+  n (8)
        //
        // ratio
        // D++ p -> p p K+ K0b	(6)
        //
        // D++ n -> p p K+ K-	(3)
        // D++ n -> p p K0 K0b	(3)
        // D++ n -> p n K+ K0b	(3)
        //
        // D+  p -> p p K+ K-	(3)
        // D+  p -> p p K0 K0b	(1)
        // D+  p -> p n K+ K0b	(3)
        //
        // D+  n -> p p K0 K-	(2)
        // D+  n -> p n K+ K-	(1)
        // D+  n -> p n K0 K0b	(3)
        // D+  n -> n n K+ K0b	(2)
        //
        
// assert((p1->isNucleon() && p2->isResonance()) || (p2->isNucleon() && p1->isResonance()));
                
        G4double sigma = 0.;
        const G4int iso = ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        const G4double ener = 0.001*KinematicsUtils::totalEnergyInCM(p1, p2); // GeV
        
        if(ener <= 2.872)
            return 0.;
        
        if(iso == 0)// D+  n
            sigma = 8* 14. * 5./19. * 0.3 *std::pow(1.-2.872*2.872/(ener*ener),3.)*std::pow(2.872*2.872/(ener*ener),0.8);
        else if (ParticleTable::getIsospin(p1->getType()) == ParticleTable::getIsospin(p2->getType()))// D+  p
            sigma = 7* 14. * 5./19. * 0.3 *std::pow(1.-2.872*2.872/(ener*ener),3.)*std::pow(2.872*2.872/(ener*ener),0.8);
        else if (std::abs(iso) == 2)// D++ n
            sigma = 9* 14. * 5./19. * 0.3 *std::pow(1.-2.872*2.872/(ener*ener),3.)*std::pow(2.872*2.872/(ener*ener),0.8);
        else // D++ p
            sigma = 6* 14. * 5./19. * 0.3 *std::pow(1.-2.872*2.872/(ener*ener),3.)*std::pow(2.872*2.872/(ener*ener),0.8);
        
        return sigma;
    }

	/// \brief piN to strange cross sections

    G4double CrossSectionsStrangeness::NpiToLK(Particle const * const p1, Particle const * const p2) {
        //
        //      Pion-Nucleon producing Lambda-Kaon cross sections
        //
        // ratio
        // p pi0 -> L K+ (1/2)
        // p pi- -> L K0 (1)
        
// assert((p1->isPion() && p2->isNucleon()) || (p2->isPion() && p1->isNucleon()));
        
        const Particle *pion;
        const Particle *nucleon;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        if(iso == 3 || iso == -3)
            return 0.;
        
        if(p1->isPion()){
            pion = p1;
            nucleon = p2;
        }
        else{
            nucleon = p1;
            pion = p2;
        }
        G4double sigma = 0.;
        
        if(pion->getType() == PiZero)
            sigma = 0.5 * p_pimToLK0(pion,nucleon);
        else
            sigma = p_pimToLK0(pion,nucleon);
        return sigma;
    }

    G4double CrossSectionsStrangeness::p_pimToLK0(Particle const * const p1, Particle const * const p2) {
        
        G4double sigma = 0.;
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(p1,p2); // GeV
        
        if(pLab < 0.911)
            return 0.;
        
        sigma = 0.3936*std::pow(pLab,-1.357)-6.052*std::exp(-std::pow(pLab-0.7154,2)/0.02026)-0.16*std::exp(-std::pow(pLab-0.9684,2)/0.001432)+0.489*std::exp(-std::pow(pLab-0.8886,2)/0.08378);
        if(sigma < 0.) return 0;
        return sigma;
    }

    G4double CrossSectionsStrangeness::NpiToSK(Particle const * const p1, Particle const * const p2) {
        //
        //      Pion-Nucleon producing Sigma-Kaon cross sections
        //
        // ratio 
        // p pi+ (5/3)    p pi0 (11/6)    p pi- (2)
        //
        // p pi+ -> S+ K+ (10)
        // p pi0 -> S+ K0 (6)*
        // p pi0 -> S0 K+ (5)
        // p pi- -> S0 K0 (6)*
        // p pi- -> S- K+ (6)
        
// assert((p1->isPion() && p2->isNucleon()) || (p2->isPion() && p1->isNucleon()));
        
        const Particle *pion;
        const Particle *nucleon;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        if(p1->isPion()){
            pion = p1;
            nucleon = p2;
        }
        else{
            nucleon = p1;
            pion = p2;
        }
        G4double sigma = 0.;
        
        if(iso == 3 || iso == -3)
            sigma = p_pipToSpKp(pion,nucleon);
        else if(pion->getType() == PiZero)
            sigma = p_pizToSzKp(pion,nucleon)+p_pimToSzKz(pion,nucleon);
        else if(iso == 1 || iso == -1)
            sigma = p_pimToSzKz(pion,nucleon)+p_pimToSmKp(pion,nucleon);
        else // should never append
            sigma = 0.;
        
        return sigma;
    }
    G4double CrossSectionsStrangeness::p_pimToSmKp(Particle const * const p1, Particle const * const p2) {
        
        G4double sigma = 0.;
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(p1,p2); // GeV
        
        if(pLab < 1.0356)
            return 0.;
        
        sigma = 4.352*std::pow(pLab-1.0356,1.006)/(std::pow(pLab+1.0356,0.0978)*std::pow(pLab,5.375));
        
        if(sigma < 0.) // should never append
            return 0;
        
        return sigma;
    }
    G4double CrossSectionsStrangeness::p_pipToSpKp(Particle const * const p1, Particle const * const p2) {
        
        G4double sigma = 0.;
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(p1,p2); // GeV
        
        if(pLab < 1.0428)
            return 0.;
        
        sigma = 0.001897*std::pow(pLab-1.0428,2.869)/(std::pow(pLab+1.0428,-16.68)*std::pow(pLab,19.1));
        
        if(sigma < 0.) // should never append
            return 0;
        
        return sigma;
    }
    G4double CrossSectionsStrangeness::p_pizToSzKp(Particle const * const p1, Particle const * const p2) {
        
        G4double sigma = 0.;
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(p1,p2); // GeV
        
        if(pLab < 1.0356)
            return 0.;
        
        sigma = 3.624*std::pow(pLab-1.0356,1.4)/std::pow(pLab,5.14);
        
        if(sigma < 0.) // should never append
            return 0;
        
        return sigma;
    }
    G4double CrossSectionsStrangeness::p_pimToSzKz(Particle const * const p1, Particle const * const p2) {
        
        G4double sigma = 0.;
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(p1,p2); // GeV
        
        if((p1->getType() == PiZero && pLab < 1.0356) || (pLab < 1.034))
            return 0.;
        
        sigma = 0.3474*std::pow(pLab-1.034,0.07678)/std::pow(pLab,1.627);
        
        if(sigma < 0.) // should never append
            return 0;
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NpiToLKpi(Particle const * const p1, Particle const * const p2) {
        //
        //      Pion-Nucleon producing Lambda-Kaon-pion cross sections
        //
        // ratio
        // p pi+ (1)    p pi0 (3/2)        p pi- (2)
        //
        // p pi0 -> L K+ pi0 (1/2)
        // all the others (1)
        //
// assert((p1->isPion() && p2->isNucleon()) || (p2->isPion() && p1->isNucleon()));
        
        G4double sigma=0.;
        const Particle *pion;
        const Particle *nucleon;
        
        if(p1->isPion()){
            pion = p1;
            nucleon = p2;
        }
        else{
            nucleon = p1;
            pion = p2;
        }
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(pion, nucleon); // GeV
        
        if(pLab < 1.147)
            return 0.;
        
        if(iso == 3 || iso == -3)
            sigma = 146.2*std::pow(pLab-1.147,1.996)/std::pow(pLab+1.147,5.921);
        else if(pion->getType() == PiZero)
            sigma = 1.5*146.2*std::pow(pLab-1.147,1.996)/std::pow(pLab+1.147,5.921);
        else
            sigma = 2*146.2*std::pow(pLab-1.147,1.996)/std::pow(pLab+1.147,5.921);
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NpiToSKpi(Particle const * const p1, Particle const * const p2) {
        //
        //      Pion-Nucleon producing Sigma-Kaon-pion cross sections
        //
        //ratio
        // p pi+ (2.25)    p pi0 (2.625)    p pi-(3)
        //
        // p pi+ -> S+ pi+ K0 (5/4)
        // p pi+ -> S+ pi0 K+ (3/4)
        // p pi+ -> S0 pi+ K+ (1/4)
        // p pi0 -> S+ pi0 K0 (1/2)
        // p pi0 -> S+ pi- K+ (1/2)
        // p pi0 -> S0 pi+ K0 (3/4)
        // p pi0 -> S0 pi0 K+ (3/8)
        // p pi0 -> S- pi+ K+ (1/2)
        // p pi- -> S+ pi- K0 (3/8)
        // p pi- -> S0 pi0 K0 (5/8)
        // p pi- -> S0 pi- K+ (5/8)
        // p pi- -> S- pi+ K0 (1)
        // p pi- -> S- pi0 K+ (3/8)
        
// assert((p1->isPion() && p2->isNucleon()) || (p2->isPion() && p1->isNucleon()));
        
        G4double sigma=0.;
        const Particle *pion;
        const Particle *nucleon;
        
        if(p1->isPion()){
            pion = p1;
            nucleon = p2;
        }
        else{
            nucleon = p1;
            pion = p2;
        }
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(pion, nucleon); // GeV
        
        if(pLab <= 1.3041)
            return 0.;
        
        if(iso == 3 || iso == -3)
            sigma = 2.25*8.139*std::pow(pLab-1.3041,2.431)/std::pow(pLab,5.298);
        else if(pion->getType() == PiZero)
            sigma = 2.625*8.139*std::pow(pLab-1.3041,2.431)/std::pow(pLab,5.298);
        else
            sigma = 3.*8.139*std::pow(pLab-1.3041,2.431)/std::pow(pLab,5.298);
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NpiToLK2pi(Particle const * const p1, Particle const * const p2) {
        //
        //      Pion-Nucleon producing Lambda-Kaon-2pion cross sections
        //
        // p pi+ (2)    p pi0 (1.75)    p pi- (2.5)
        //
        // p pi+ -> L K+ pi+ pi0 (1)
        // p pi+ -> L K0 pi+ pi+ (1)
        // p pi0 -> L K+ pi0 pi0 (1/4)
        // p pi0 -> L K+ pi+ pi- (1)
        // p pi0 -> L K0 pi+ pi0 (1/2)
        // p pi- -> L K+ pi0 pi- (1)
        // p pi- -> L K0 pi+ pi- (1)
        // p pi- -> L K0 pi0 pi0 (1/2)
        
// assert((p1->isPion() && p2->isNucleon()) || (p2->isPion() && p1->isNucleon()));
        
        G4double sigma=0.;
        const Particle *pion;
        const Particle *nucleon;
        
        if(p1->isPion()){
            pion = p1;
            nucleon = p2;
        }
        else{
            nucleon = p1;
            pion = p2;
        }
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(pion, nucleon); // GeV
        
        if(pLab <= 1.4162)
            return 0.;
        
        if(iso == 3 || iso == -3)
            sigma = 2*18.77*std::pow(pLab-1.4162,4.597)/std::pow(pLab,6.877);
        else if(pion->getType() == PiZero)
            sigma = 1.75*18.77*std::pow(pLab-1.4162,4.597)/std::pow(pLab,6.877);
        else
            sigma = 2.5*18.77*std::pow(pLab-1.4162,4.597)/std::pow(pLab,6.877);
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NpiToSK2pi(Particle const * const p1, Particle const * const p2) {
        //
        //      Pion-Nucleon producing Lambda-Kaon-2pion cross sections
        //
        // ratio
        // p pi+ (3.25)    p pi0 (3.5)    p pi- (3.75)
        //
        // p pi+ -> S+ K+ pi+ pi- (1)
        // p pi+ -> S+ K+ pi0 pi0 (1/4)
        // p pi+ -> S0 K+ pi+ pi0 (1/2)
        // p pi+ -> S- K+ pi+ pi+ (1/4)
        // p pi+ -> S+ K0 pi+ pi0 (1)
        // p pi+ -> S0 K0 pi+ pi+ (1/4)
        //
        // p pi0 -> S+ K+ pi0 pi- (1/2)
        // p pi0 -> S0 K+ pi+ pi- (1/2)
        // p pi0 -> S0 K+ pi0 pi0 (1/4)
        // p pi0 -> S- K+ pi+ pi0 (1/4)
        // p pi0 -> S+ K0 pi+ pi- (1)
        // p pi0 -> S+ K0 pi0 pi0 (1/4)
        // p pi0 -> S0 K0 pi+ pi0 (1/4)
        // p pi0 -> S- K0 pi+ pi+ (1/2)
        //
        // p pi- -> S+ K+ pi- pi- (1/4)
        // p pi- -> S0 K+ pi0 pi- (1/2)
        // p pi- -> S- K+ pi+ pi- (1/4)
        // p pi- -> S- K+ pi0 pi0 (1/4)
        // p pi- -> S+ K0 pi0 pi- (1/2)
        // p pi- -> S0 K0 pi+ pi- (1)
        // p pi- -> S0 K0 pi0 pi0 (1/2)
        // p pi- -> S- K0 pi+ pi0 (1/2)
        
// assert((p1->isPion() && p2->isNucleon()) || (p2->isPion() && p1->isNucleon()));
        
        G4double sigma=0.;
        const Particle *pion;
        const Particle *nucleon;
        
        if(p1->isPion()){
            pion = p1;
            nucleon = p2;
        }
        else{
            nucleon = p1;
            pion = p2;
        }
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(pion, nucleon); // GeV
        
        if(pLab <= 1.5851)
            return 0.;
        
        if(iso == 3 || iso == -3)
            sigma = 3.25*137.6*std::pow(pLab-1.5851,5.856)/std::pow(pLab,9.295);
        else if(pion->getType() == PiZero)
            sigma = 3.5*137.6*std::pow(pLab-1.5851,5.856)/std::pow(pLab,9.295);
        else
            sigma = 3.75*137.6*std::pow(pLab-1.5851,5.856)/std::pow(pLab,9.295);
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NpiToNKKb(Particle const * const p1, Particle const * const p2) {
        //
        //      Pion-Nucleon producing Nucleon-Kaon-antiKaon cross sections
        //
        // ratio
        // p pi+ (1/2)    p pi0 (3/2)    p pi- (5/2)    
        //
        // p pi+ -> p K+ K0b (1/2)
        // p pi0 -> p K0 K0b (1/4)
        // p pi0 -> p K+ K- (1/4)
        // p pi0 -> n K+ K0b (1)
        // p pi- -> p K0 K- (1/2)
        // p pi- -> n K+ K- (1)
        // p pi- -> n K0 K0b (1)
        
// assert((p1->isPion() && p2->isNucleon()) || (p2->isPion() && p1->isNucleon()));
        
        const Particle *particle1;
        const Particle *particle2;
        
        if(p1->isPion()){
            particle1 = p1;
            particle2 = p2;
        }
        else{
            particle1 = p2;
            particle2 = p1;
        }
        
        G4double sigma = 0.;
        
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(particle1,particle2); // GeV
        
        if(particle1->getType() == PiZero){
            if(pLab < 1.5066) return 0.;
            else if(pLab < 30.) sigma = 3./2.*2.996*std::pow((pLab - 1.5066),1.929)/std::pow(pLab,3.582);
            else return 0.;
        }
        else if((particle1->getType() == PiPlus && particle2->getType() == Neutron) || (particle1->getType() == PiMinus && particle2->getType() == Proton)){
            if(pLab < 1.5066) return 0.;
            else if(pLab < 30.) sigma = 5./2.*2.996*std::pow((pLab - 1.5066),1.929)/std::pow(pLab,3.582);
            else return 0.;
        }
        else{
            if(pLab < 1.5066) return 0.;
            else if(pLab < 30.) sigma = 1./2.*2.996*std::pow((pLab - 1.5066),1.929)/std::pow(pLab,3.582);
            else return 0.;
        }
        return sigma;
    }

    G4double CrossSectionsStrangeness::NpiToMissingStrangeness(Particle const * const p1, Particle const * const p2) {
        //
        //      Pion-Nucleon missing strangeness production cross sections
        //
// assert((p1->isPion() && p2->isNucleon()) || (p2->isPion() && p1->isNucleon()));
        
        const Particle *pion;
        const Particle *nucleon;
        
        if(p1->isPion()){
            pion = p1;
            nucleon = p2;
        }
        else{
            pion = p2;
            nucleon = p1;
        }
        
        G4double sigma = 0.;
        
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(pion,nucleon); // GeV
        if(pLab < 2.2) return 0.;
        
        if(pion->getType() == PiZero){
            if(pLab < 30.) sigma = 4.4755*std::pow((pLab - 2.2),1.927)/std::pow(pLab,1.89343);
            else return 0.;
        }
        else if((pion->getType() == PiPlus && nucleon->getType() == Neutron) || (pion->getType() == PiMinus && nucleon->getType() == Proton)){
            if(pLab < 30.) sigma = 5.1*std::pow((pLab - 2.2),1.854)/std::pow(pLab,1.904);
            else return 0.;
        }
        else{
            if(pLab < 30.) sigma = 3.851*std::pow((pLab - 2.2),2)/std::pow(pLab,1.88286);
            else return 0.;
        }
        return sigma;
    }

	/// \brief NY cross sections

    G4double CrossSectionsStrangeness::NLToNS(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Lambda producing Nucleon-Sigma cross sections
        //
        // ratio
        // p L -> p S0 (1/2)
        // p L -> n S+ (1)
        
        
// assert((p1->isLambda() && p2->isNucleon()) || (p2->isLambda() && p1->isNucleon()));
        
        G4double sigma = 0.;
        
        const Particle *particle1;
        const Particle *particle2;
        
        if(p1->isLambda()){
            particle1 = p1;
            particle2 = p2;
        }
        else{
            particle1 = p2;
            particle2 = p1;
        }
        
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(particle1,particle2);
        
        if(pLab < 0.664)
            return 0.;
        
        sigma = 3 * 8.74*std::pow((pLab-0.664),0.438)/std::pow(pLab,2.717); // 3 * L p -> S0 p
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NSToNL(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Lambda producing Nucleon-Sigma cross sections
        //
        // ratio
        // p S0 -> p L (1/2)
        // p S- -> n L (1)
        
// assert((p1->isSigma() && p2->isNucleon()) || (p2->isSigma() && p1->isNucleon()));
        
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        if(iso == 3 || iso == -3)
            return 0.;
        
        G4double sigma;
        const Particle *particle1;
        const Particle *particle2;
        
        if(p1->isSigma()){
            particle1 = p1;
            particle2 = p2;
        }
        else{
            particle1 = p2;
            particle2 = p1;
        }
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(particle1,particle2); //  GeV
        
        if(particle1->getType() == SigmaZero){
            if(pLab < 0.1) return 100.; // cut-max
            sigma = 8.23*std::pow(pLab,-1.087);
        }
        else{
            if(pLab < 0.1) return 200.; // cut-max
            sigma = 16.46*std::pow(pLab,-1.087);
        }
        return sigma;
    }
    
    G4double CrossSectionsStrangeness::NSToNS(Particle const * const p1, Particle const * const p2) {
        
// assert((p1->isSigma() && p2->isNucleon()) || (p2->isSigma() && p1->isNucleon()));
        
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        if(iso == 3 || iso == -3)
            return 0.; // only quasi-elastic here
        
        const Particle *particle1;
        const Particle *particle2;
        
        if(p1->isSigma()){
            particle1 = p1;
            particle2 = p2;
        }
        else{
            particle1 = p2;
            particle2 = p1;
        }
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(particle1,particle2); //  GeV
        
        if(particle2->getType() == Neutron && pLab < 0.162) return 0.;
        else if(pLab < 0.1035) return 200.; // cut-max
        
        return 13.79*std::pow(pLab,-1.181);
    }

	/// \brief NK cross sections

    G4double CrossSectionsStrangeness::NKToNK(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Kaon quasi-elastic cross sections
        //
// assert((p1->isNucleon() && p2->isKaon()) || (p2->isNucleon() && p1->isKaon()));
        
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        if(iso != 0)
            return 0.;
        
        const Particle *particle1;
        const Particle *particle2;
        
        if(p1->isKaon()){
            particle1 = p1;
            particle2 = p2;
        }
        else{
            particle1 = p2;
            particle2 = p1;
        }
        
        G4double sigma = 0.;
        G4double pLab = 0.001 * KinematicsUtils::momentumInLab(particle1,particle2); // GeV
        
        if(particle1->getType() == Proton)
            pLab += 2*0.0774;
        
        if(pLab <= 0.0774)
            return 0.;
        
        sigma = 12.84*std::pow((pLab-0.0774),18.19)/std::pow((pLab),20.41);
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NKToNKpi(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Kaon producing Nucleon-Kaon-pion cross sections
        //
        // Ratio determined by meson symmetry using only "resonante" diagram (with Delta or K*)
        //
        // ratio: K+ p (5)    K0 p (5.545)
        //
        // K+ p -> p K+ pi0        1.2
        // K+ p -> p K0 pi+        3
        // K+ p -> n K+ pi+        0.8
        // K0 p -> p K+ pi-        1
        // K0 p -> p K0 pi0        0.845
        // K0 p -> n K+ pi0        1.47
        // K0 p -> n K0 pi+        2.23
// assert((p1->isNucleon() && p2->isKaon()) || (p2->isNucleon() && p1->isKaon()));
        
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const Particle *particle1;
        const Particle *particle2;
        
        if(p1->isKaon()){
            particle1 = p1;
            particle2 = p2;
        }
        else{
            particle1 = p2;
            particle2 = p1;
        }
        
        G4double sigma = 0.;
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(particle1,particle2); // GeV
        
        if(pLab <= 0.53)
            return 0.;
        
        if(iso == 0)
            sigma = 5.55*116.8*std::pow(pLab-0.53,6.874)/std::pow(pLab,10.11);
        else
            sigma = 5.*116.8*std::pow(pLab-0.53,6.874)/std::pow(pLab,10.11);;
        return sigma;
    }

    G4double CrossSectionsStrangeness::NKToNK2pi(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-Kaon producing Nucleon-Kaon-2pion cross sections
        //
        // p K+ (2.875) p K0 (3.125)
        //
        // p K+ -> p K+ pi+ pi- (1)
        // p K+ -> p K+ pi0 pi0 (1/8)
        // p K+ -> p K0 pi+ pi0 (1)
        // p K+ -> n K+ pi+ pi0 (1/2)
        // p K+ -> n K0 pi+ pi+ (1/4)
        // p K0 -> p K+ pi0 pi- (1)
        // p K0 -> p K0 pi+ pi- (1)
        // p K0 -> p K0 pi0 pi0 (1/8)
        // p K0 -> n K+ pi+ pi- (1/4)
        // p K0 -> n K+ pi0 pi0 (1/4)
        // p K0 -> n K0 pi+ pi0 (1/2)
        
// assert((p1->isNucleon() && p2->isKaon()) || (p2->isNucleon() && p1->isKaon()));
        
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const Particle *particle1;
        const Particle *particle2;
        
        if(p1->isKaon()){
            particle1 = p1;
            particle2 = p2;
        }
        else{
            particle1 = p2;
            particle2 = p1;
        }
        
        G4double sigma = 0.;
        const G4double pLab = 0.001 * KinematicsUtils::momentumInLab(particle1,particle2); // GeV
        
        if(pLab < 0.812)
            sigma = 0.;
        else if(pLab < 1.744)
            sigma = 26.41*std::pow(pLab-0.812,7.138)/std::pow(pLab,5.337);
        else if(pLab < 3.728)
            sigma = 1572.*std::pow(pLab-0.812,9.069)/std::pow(pLab,12.44);
        else
            sigma = 60.23*std::pow(pLab-0.812,5.084)/std::pow(pLab,6.72);
        
        if(iso == 0)
            sigma *= 3.125;
        else
            sigma *= 2.875;
        
        return sigma;
    }

	/// \brief NKb cross sections

    G4double CrossSectionsStrangeness::NKbToNKb(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-antiKaon quasi-elastic cross sections
        //
// assert((p1->isNucleon() && p2->isAntiKaon()) || (p1->isAntiKaon() && p2->isNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const Particle *antikaon;
        const Particle *nucleon;
        
        if (p1->isAntiKaon()) {
            antikaon = p1;
            nucleon = p2;
        }
        else {
            antikaon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antikaon, nucleon); // GeV
        
        if(iso != 0) // K0b p and K- n -> forbidden: quasi-elastic diffusion only
            return 0;
        else if(nucleon->getType() == Proton){ // K- p -> K0b n
            if(pLab < 0.08921)
                return 0.;
            else if(pLab < 0.2)
                sigma = 0.4977*std::pow(pLab - 0.08921,0.5581)/std::pow(pLab,2.704);
            else if(pLab < 0.73)
                sigma = 2.*std::pow(pLab,-1.2) + 6.493*std::exp(-0.5*std::pow((pLab-0.3962)/0.02,2));
            else if(pLab < 1.38)
                sigma = 2.3*std::pow(pLab,-0.9) + 1.1*std::exp(-0.5*std::pow((pLab-0.82)/0.04,2)) + 5.*std::exp(-0.5*std::pow((pLab-1.04)/0.1,2));
            else if(pLab < 30.)
                sigma = 2.5*std::pow(pLab,-1.68) + 0.7*std::exp(-0.5*std::pow((pLab-1.6)/0.2,2)) + 0.2*std::exp(-0.5*std::pow((pLab-2.3)/0.2,2));
            else sigma = 0.;
        }
        else{ // K0b n -> K- p (same as K- p but without threshold)
            if(pLab < 0.1)
                sigma = 30.;
            else if(pLab < 0.73)
                sigma = 2.*std::pow(pLab,-1.2) + 6.493*std::exp(-0.5*std::pow((pLab-0.3962)/0.02,2));
            else if(pLab < 1.38)
                sigma = 2.3*std::pow(pLab,-0.9) + 1.1*std::exp(-0.5*std::pow((pLab-0.82)/0.04,2)) + 5.*std::exp(-0.5*std::pow((pLab-1.04)/0.1,2));
            else if(pLab < 30.)
                sigma = 2.5*std::pow(pLab,-1.68) + 0.7*std::exp(-0.5*std::pow((pLab-1.6)/0.2,2)) + 0.2*std::exp(-0.5*std::pow((pLab-2.3)/0.2,2));
            else sigma = 0.;
        }
        return sigma;
    }

    G4double CrossSectionsStrangeness::NKbToSpi(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-antiKaon producing Sigma-pion cross sections
        //
        // ratio
        // p K0b (4/3)    p K- (13/6)
        //
        // p K0b -> S+ pi0 (2/3)
        // p K0b -> S0 pi+ (2/3)
        // p K-  -> S+ pi- (1)
        // p K-  -> S0 pi0 (1/2)
        // p K-  -> S- pi+ (2/3)
        
// assert((p1->isNucleon() && p2->isAntiKaon()) || (p1->isAntiKaon() && p2->isNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const Particle *antikaon;
        const Particle *nucleon;
        
        if (p1->isAntiKaon()) {
            antikaon = p1;
            nucleon = p2;
        }
        else {
            antikaon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antikaon, nucleon); // GeV
        
        if(iso == 0){
            if(pLab < 0.1)
                return 152.0; // 70.166*13/6
            else
                sigma = 13./6.*(1.4*std::pow(pLab,-1.7)+1.88*std::exp(-std::pow(pLab-0.747,2)/0.005)+8*std::exp(-std::pow(pLab-0.4,2)/0.002)+0.8*std::exp(-std::pow(pLab-1.07,2)/0.01));
        }
        else{
            if(pLab < 0.1)
                return 93.555; // 70.166*4/3
            else
                sigma = 4./3.*(1.4*std::pow(pLab,-1.7)+1.88*std::exp(-std::pow(pLab-0.747,2)/0.005)+8*std::exp(-std::pow(pLab-0.4,2)/0.002)+0.8*std::exp(-std::pow(pLab-1.07,2)/0.01));
                //sigma = 4./3.*(1.4*std::pow(pLab,-1.7)+1.88*std::exp(-std::pow(pLab-0.747,2)/0.005)+0.8*std::exp(-std::pow(pLab-1.07,2)/0.01));
        }
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NKbToLpi(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-antiKaon producing Lambda-pion cross sections
        //
        // ratio
        // p K0b (1)    p K- (1/2)
        //
        // p K- -> L pi0 (1/2)
        // p K0b -> L pi+ (1)
// assert((p1->isNucleon() && p2->isAntiKaon()) || (p1->isAntiKaon() && p2->isNucleon()));
        
        G4double sigma = 0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const Particle *antikaon;
        const Particle *nucleon;
        
        if (p1->isAntiKaon()) {
            antikaon = p1;
            nucleon = p2;
        }
        else {
            antikaon = p2;
            nucleon = p1;
        }
        if(iso == 0)
            sigma = p_kmToL_pz(antikaon,nucleon);
        else
            sigma = 2*p_kmToL_pz(antikaon,nucleon);
        
        return sigma;
    }
    G4double CrossSectionsStrangeness::p_kmToL_pz(Particle const * const p1, Particle const * const p2) {
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(p1, p2); // GeV
        G4double sigma = 0.;
        if(pLab < 0.086636)
            sigma = 40.24;
        else if(pLab < 0.5)
            sigma = 0.97*std::pow(pLab,-1.523);
        else if(pLab < 2.)
            sigma = 1.23*std::pow(pLab,-1.467)+0.872*std::exp(-std::pow(pLab-0.749,2)/0.0045)+2.337*std::exp(-std::pow(pLab-0.957,2)/0.017)+0.476*std::exp(-std::pow(pLab-1.434,2)/0.136);
        else if(pLab < 30.)
            sigma = 3.*std::pow(pLab,-2.57);
        else
            sigma = 0.;
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NKbToS2pi(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-antiKaon producing Sigma-2pion cross sections
        //
        // ratio
        // p K0b (29/12)    p K- (59/24)
        //
        // p K0b -> S+ pi+ pi- (2/3)
        // p K0b -> S+ pi0 pi0 (1/4)
        // p K0b -> S0 pi+ pi0 (5/6)
        // p K0b -> S- pi+ pi+ (2/3)
        // p K-  -> S+ pi0 pi- (1)
        // p K-  -> S0 pi+ pi- (2/3)
        // p K-  -> S0 pi0 pi0 (1/8)
        // p K-  -> S- pi+ pi0 (2/3)
        
// assert((p1->isNucleon() && p2->isAntiKaon()) || (p1->isAntiKaon() && p2->isNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const Particle *antikaon;
        const Particle *nucleon;
        
        if (p1->isAntiKaon()) {
            antikaon = p1;
            nucleon = p2;
        }
        else {
            antikaon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antikaon, nucleon); // GeV
        
        if(pLab < 0.260)
            return 0.;
        
        if(iso == 0)
            sigma = 29./12.*3./2.*(49.96*std::pow(pLab-0.260,6.398)/std::pow(pLab+0.260,9.732)+0.1451*std::exp(-std::pow(pLab-0.4031,2)/0.00115));
        else
            sigma = 54./24.*3./2.*(49.96*std::pow(pLab-0.260,6.398)/std::pow(pLab+0.260,9.732)+0.1451*std::exp(-std::pow(pLab-0.4031,2)/0.00115));
        
        /*
        if(iso == 0)
            sigma = 29./12.*(85.46*std::pow(pLab-0.226,8.118)/std::pow(pLab+0.226,11.69)+0.1451*std::exp(-std::pow(pLab-0.4031,2)/0.00115));
        else
            sigma = 54./24.*(85.46*std::pow(pLab-0.226,8.118)/std::pow(pLab+0.226,11.69)+0.1451*std::exp(-std::pow(pLab-0.4031,2)/0.00115));*/
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NKbToL2pi(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-antiKaon producing Lambda-2pion cross sections
        //
        // ratio
        // p K0b -> L pi+ pi0 (1)
        // p K- -> L pi+ pi- (1)
        // p K- -> L pi0 pi0 (1/4)
        
// assert((p1->isNucleon() && p2->isAntiKaon()) || (p1->isAntiKaon() && p2->isNucleon()));
        
        G4double sigma = 0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const Particle *antikaon;
        const Particle *nucleon;
        
        if (p1->isAntiKaon()) {
            antikaon = p1;
            nucleon = p2;
        }
        else {
            antikaon = p2;
            nucleon = p1;
        }
        
        if(iso == 0)
            sigma = 1.25*p_kmToL_pp_pm(antikaon,nucleon);
        else
            sigma = p_kmToL_pp_pm(antikaon,nucleon);
        
        return sigma;
    }
    G4double CrossSectionsStrangeness::p_kmToL_pp_pm(Particle const * const p1, Particle const * const p2) {
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(p1, p2); // GeV
        G4double sigma = 0.;
        if(pLab < 0.97)
            sigma = 6364.*std::pow(pLab,6.07)/std::pow(pLab+1.,10.58)+2.158*std::exp(-std::pow((pLab-0.395)/.01984,2)/2.);
        else if(pLab < 30)
            sigma = 46.3*std::pow(pLab,0.62)/std::pow(pLab+1.,3.565);
        else
            sigma = 0.;
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NKbToNKbpi(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-antiKaon producing Nucleon-antiKaon-pion cross sections
        //
        // ratio
        // p K0b (2)    p K- (4)
        //
        // p K0b -> p K0b pi0 (1/2)
        // p K0b -> p K- pi+ (1)*
        // p K0b -> n K0b pi+ (1/2)*
        // p K- -> p K- pi0 (1/2)*
        // p K- -> p K0b pi- (2/3)*
        // p K- -> n K- pi+ (3/4)*
        // p K- -> n K0b pi0 (2)
        
        // 
// assert((p1->isNucleon() && p2->isAntiKaon()) || (p1->isAntiKaon() && p2->isNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const Particle *antikaon;
        const Particle *nucleon;
        
        if (p1->isAntiKaon()) {
            antikaon = p1;
            nucleon = p2;
        }
        else {
            antikaon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antikaon, nucleon); // GeV
        
        if(pLab < 0.526)
            return 0.;
        
        if(iso == 0)
            sigma = 2. * 101.3*std::pow(pLab-0.526,5.846)/std::pow(pLab,8.343);
        else
            sigma = 4. * 101.3*std::pow(pLab-0.526,5.846)/std::pow(pLab,8.343);
        
        return sigma;
    }

    G4double CrossSectionsStrangeness::NKbToNKb2pi(Particle const * const p1, Particle const * const p2) {
        //
        //      Nucleon-antiKaon producing Nucleon-antiKaon-2pion cross sections
        //
        // ratio
        // p K0b (4.25)    p K- (4.75)
        //
        // p K0b -> p K0b pi+ pi- (1)
        // p K0b -> p K0b pi0 pi0 (1/4)
        // p K0b -> p K-  pi+ pi0 (1)
        // p K0b -> n K0b pi+ pi0 (1)
        // p K0b -> n K-  pi+ pi+ (1)
        // p K-  -> p K0b pi0 pi- (1)
        // p K-  -> p K-  pi+ pi- (1)
        // p K-  -> p K-  pi0 pi0 (1/4)
        // p K-  -> n K0b pi+ pi- (1)
        // p K-  -> n K0b pi0 pi0 (1/2)
        // p K-  -> n K-  pi+ pi0 (1)
        
// assert((p1->isNucleon() && p2->isAntiKaon()) || (p1->isAntiKaon() && p2->isNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        
        const Particle *antikaon;
        const Particle *nucleon;
        
        if (p1->isAntiKaon()) {
            antikaon = p1;
            nucleon = p2;
        }
        else {
            antikaon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antikaon, nucleon); // GeV
        
        if(pLab < 0.85)
            return 0.;
        
        if(iso == 0)
            sigma = 4.75 * 26.8*std::pow(pLab-0.85,4.9)/std::pow(pLab,6.34);
        else
            sigma = 4.25 * 26.8*std::pow(pLab-0.85,4.9)/std::pow(pLab,6.34);
        
        return sigma;
    }


} // namespace G4INCL

