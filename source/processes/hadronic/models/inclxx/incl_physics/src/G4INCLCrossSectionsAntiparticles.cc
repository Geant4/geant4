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

/** \file G4INCLCrossSectionsAntiparticles.cc
 * \brief Multipion, mesonic Resonances, strange cross sections and antinucleon as projectile
 *
 * \date 31st March 2023
 * \author Demid Zharenov
 */

#include "G4INCLCrossSectionsAntiparticles.hh"
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
    
    const G4int CrossSectionsAntiparticles::nMaxPiNN = 4;
    const G4int CrossSectionsAntiparticles::nMaxPiPiN = 4;

    CrossSectionsAntiparticles::CrossSectionsAntiparticles() :
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

    G4double CrossSectionsAntiparticles::total(Particle const * const p1, Particle const * const p2) {
        G4double inelastic;
        if ((p1->isNucleon() && p2->isAntiNucleon()) || (p1->isAntiNucleon() && p2->isNucleon()))
            inelastic = NNbarCEX(p1, p2) + NNbarToNNbarpi(p1, p2) + NNbarToNNbar2pi(p1, p2) + NNbarToNNbar3pi(p1, p2) + NNbarToAnnihilation(p1, p2) + NNbarToLLbar(p1, p2);   
        else if(p1->isNucleon() && p2->isNucleon()) {
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
            inelastic = CrossSectionsStrangeness::NLToNS(p1,p2);
        } else if((p1->isNucleon() && p2->isSigma()) ||
                  (p1->isSigma() && p2->isNucleon())) {
            inelastic = CrossSectionsStrangeness::NSToNL(p1,p2) + CrossSectionsStrangeness::NSToNS(p1,p2);
        } else if((p1->isNucleon() && p2->isKaon()) ||
                  (p1->isKaon() && p2->isNucleon())) {
            inelastic = CrossSectionsStrangeness::NKToNK(p1,p2) + CrossSectionsStrangeness::NKToNKpi(p1,p2) + CrossSectionsStrangeness::NKToNK2pi(p1,p2);
        } else if((p1->isNucleon() && p2->isAntiKaon()) ||
                  (p1->isAntiKaon() && p2->isNucleon())) {
            inelastic = CrossSectionsStrangeness::NKbToLpi(p1,p2) 
            + CrossSectionsStrangeness::NKbToSpi(p1,p2) + CrossSectionsStrangeness::NKbToL2pi(p1,p2) 
            + CrossSectionsStrangeness::NKbToS2pi(p1,p2) + CrossSectionsStrangeness::NKbToNKb(p1,p2) 
            + CrossSectionsStrangeness::NKbToNKbpi(p1,p2) + CrossSectionsStrangeness::NKbToNKb2pi(p1,p2);
        } else {
            inelastic = 0.;
        }
        return inelastic + elastic(p1, p2);
    }    

    // without NNbar!
    G4double CrossSectionsAntiparticles::elastic(Particle const * const p1, Particle const * const p2) {
        if ((p1->isNucleon() && p2->isAntiNucleon()) || (p1->isAntiNucleon() && p2->isNucleon()))
            return NNbarElastic(p1, p2);
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
            return CrossSectionsStrangeness::NYelastic(p1, p2);
        }
        else if ((p1->isNucleon() && p2->isKaon()) || (p2->isNucleon() && p1->isKaon())){
            return CrossSectionsStrangeness::NKelastic(p1, p2);
        }
        else if ((p1->isNucleon() && p2->isAntiKaon()) || (p2->isNucleon() && p1->isAntiKaon())){
            return CrossSectionsStrangeness::NKbelastic(p1, p2);
        }
        else {
            return 0.0;
        }
    }

    G4double CrossSectionsAntiparticles::NNbarCEX(Particle const * const p1, Particle const * const p2) {
        //brief ppbar
        // p pbar -> n nbar (BFMM 204)
        //
        //brief nnbar
        // n nbar -> p pbar (same as BFMM 204, but no threshold)
        //

// assert((p1->isAntiNucleon() && p2->isNucleon()) || (p1->isNucleon() && p2->isAntiNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        // iso == 2 || iso == -2 (n pbar or p nbar)

        const std::vector<G4double> BFMM204 = {7.549, -0.041, -2.959, -6.835, 1.629, 0.114};
        //{6.875, 0.590, -0.003, -6.629, 1.532, 0.114}
        //const G4double Eth_PPbar_NNbar = 0.114;
        const std::vector<G4double> BFMM204nn = {7.549, -0.041, -2.959, -6.835, 1.629};
        //const G4double Eth_NNbar_PPbar = 0.0;

        const Particle *antinucleon;
        const Particle *nucleon;
        
        if (p1->isAntiNucleon()) {
            antinucleon = p1;
            nucleon = p2;
        }
        else {
            antinucleon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antinucleon, nucleon); // GeV
        
        if(iso == 2 || iso == -2){ //npbar or pnbar
            sigma = 0.0;
            return sigma;
        }
        else{ // ppbar or nnbar
            if(p1->getType()==antiProton || p1->getType()==Proton)
                sigma = KinematicsUtils::compute_xs(BFMM204, pLab); // ppbar case
            else
                sigma = KinematicsUtils::compute_xs(BFMM204nn, pLab); // nnbar case
            return sigma;
        }
    }

    G4double CrossSectionsAntiparticles::NNbarElastic(Particle const * const p1, Particle const * const p2) {
        //brief ppbar
        // p pbar -> p pbar (BFMM 2)
        //
        //brief npbar
        // n pbar -> n pbar (BFMM 472)
        //
        //brief nnbar
        // n nbar -> n nbar (same as BFMM 2)
        //
        //brief pnbar 
        // p nbar -> p nbar (same as BFMM 472)
        //

// assert((p1->isAntiNucleon() && p2->isNucleon()) || (p1->isNucleon() && p2->isAntiNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        // iso == 2 || iso == -2 (n pbar or p nbar)

        const std::vector<G4double> BFMM2 = {110.496, -65.605, -0.198, -34.813, 4.317};
        //elastic ppbar;
        const std::vector<G4double> BFMM472 = {14.625, 23.413, -0.288, -9.002, 1.084};
        //elastic pnbar;

        const Particle *antinucleon;
        const Particle *nucleon;
        
        if (p1->isAntiNucleon()) {
            antinucleon = p1;
            nucleon = p2;
        }
        else {
            antinucleon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antinucleon, nucleon); // GeV
        
        if(iso == 2 || iso == -2){ //npbar or pnbar
            sigma = KinematicsUtils::compute_xs(BFMM472, pLab);
            return sigma;
        }
        else{ // ppbar or nnbar
            if(p1->getType()==antiProton || p1->getType()==Proton)
                sigma = KinematicsUtils::compute_xs(BFMM2, pLab); // ppbar case
            else
                sigma = KinematicsUtils::compute_xs(BFMM2, pLab); // nnbar case
            return sigma;
        }
    }

    G4double CrossSectionsAntiparticles::NNbarToLLbar(Particle const * const p1, Particle const * const p2) {
        // this channel includes all states with lambdas, sigmas and xis and their antiparticles

        //brief ppbar
        // p pbar -> l lbar (BFMM 121)
        // ppbar -> l lbar pi0 (BFMM 113)
        // ppbar -> splus pim lbar || sminusbar pim l (BFMM 136)
        // ppbar -> sminus pip lbar || splusbar l pip (BFMM 146)
        // ppbar -> sp spbar (BFMM 139)
        // ppbar -> sm smbar (BFMM 149) 
        // ppbar -> szero szerobar (BFMM 144)
        // ppbar -> ximinus ximinusbar (BFMM 101)
        // ppbar -> szero lbar || szerobar l (BFMM 143)
        //
        //
        //brief npbar
        // n pbar -> l lbar pi- (BFMM 487)
        // n pbar -> l sbarplus || lbar sminus (BFMM 488)
        //
        //
        //brief nnbar
        // all same as for ppbar
        //
        //
        //brief pnbar
        // p nbar -> l lbar pi+ (same as BFMM 487)
        // p nbar -> l sbarminus || lbar splus (same as BFMM 488)
        //

        const std::vector<G4double> BFMM121 = {2.379, -2.738, -1.260, -1.915, 0.430, 1.437};
        //const G4double Eth_PPbar_LLbar = 1.437;
        const std::vector<G4double> BFMM113 = {-0.105, 0.000, -5.099, 0.188, -0.050, 1.820};
        //const G4double Eth_PPbar_LLbar_pi0 = 1.820;
        const std::vector<G4double> BFMM139 = {0.142, -0.291, -1.702, -0.058, 0.001, 1.851};
        //const G4double Eth_PPbar_SpSpbar = 1.851;
        const std::vector<G4double> BFMM149 = {1.855, -2.238, -1.002, -1.279, 0.252, 1.896};
        //const G4double Eth_PPbar_SmSmbar = 1.896;
        const std::vector<G4double> BFMM136 = {1.749, -2.506, -1.222, -1.262, 0.274, 2.042};
        //const G4double Eth_PPbar_SpLbar_pim = 2.042;
        const std::vector<G4double> BFMM146 = {1.037, -1.437, -1.155, -0.709, 0.138, 2.065};
        //const G4double Eth_PPbar_SmLbar_pip = 2.065;
        const std::vector<G4double> BFMM143 = {0.652, -1.006, -1.805, -0.537, 0.121, 1.653};
        //const G4double Eth_PPbar_Szero_Lbar = 1.653;

// assert((p1->isAntiNucleon() && p2->isNucleon()) || (p1->isNucleon() && p2->isAntiNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        // iso == 2 || iso == -2 (n pbar or p nbar)

        const Particle *antinucleon;
        const Particle *nucleon;
        
        if (p1->isAntiNucleon()) {
            antinucleon = p1;
            nucleon = p2;
        }
        else {
            antinucleon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antinucleon, nucleon); // GeV

        //fixed due to limited data
        G4double BFMM144; 
        if(pLab > 1.868) BFMM144 = 0.008; //sigmazero sigmazerobar
        else BFMM144 = 0.0;
        G4double BFMM101; 
        if(pLab > 1.868) BFMM101 = 0.002; //xizero xizerobar
        else BFMM101 = 0.0;

        // npbar cross sections (fixed due to limited data)
        G4double BFMM487;
        if(pLab > 2.1) BFMM487 = 0.048; //llbar piminus
        else BFMM487 = 0.0;
        G4double BFMM488;
        if(pLab > 2.0) BFMM488 = 0.139; //lsigmaminus +cc
        else BFMM488 = 0.0;
        
        if(iso == 2 || iso == -2){ //npbar or pnbar
            sigma = BFMM487 + BFMM488;
            return sigma;
        }
        else{ // ppbar or nnbar
            sigma = KinematicsUtils::compute_xs(BFMM113, pLab) 
            +KinematicsUtils::compute_xs(BFMM139, pLab) +KinematicsUtils::compute_xs(BFMM136, pLab)
            +KinematicsUtils::compute_xs(BFMM146, pLab)+KinematicsUtils::compute_xs(BFMM143, pLab) 
            +KinematicsUtils::compute_xs(BFMM121, pLab)+KinematicsUtils::compute_xs(BFMM149, pLab)
            +BFMM144 +BFMM101; // nnbar case totally same as ppbar
            return sigma;
        }
    }

    G4double CrossSectionsAntiparticles::NNbarToNNbarpi(Particle const * const p1, Particle const * const p2) {
        //brief ppbar
        // p pbar -> p pbar pi0 (BFMM 185)
        // p pbar -> p nbar pi- (BFMM 188)
        // p pbar -> n pbar pi+ (BFMM 199)
        // p pbar -> n nbar pi0 (no data)
        //
        //brief npbar
        // n pbar -> p pbar pi- (BFMM 491)
        // n pbar -> p nbar pion (impossible)
        // n pbar -> n pbar pi0 (BFMM 495)
        // n pbar -> n nbar pi- (same as BFMM 188)
        //
        //brief nnbar
        // n nbar -> n nbar pi0 (same as BFMM 185)
        // n nbar -> p nbar pi- (same as BFMM 188)
        // n nbar -> n pbar pi+ (same as BFMM 199)
        // n nbar -> p pbar pi0 (no data)
        //
        //brief pnbar 
        // p nbar -> p pbar pi+ (same as BFMM 491)
        // p nbar -> n pbar pion (impossible)
        // p nbar -> p nbar pi0 (BFMM 495)
        // p nbar -> n nbar pi- (same as BFMM 188)
        //
        //
        // BFMM 188,199 are very close in value, 491 is larger

// assert((p1->isAntiNucleon() && p2->isNucleon()) || (p1->isNucleon() && p2->isAntiNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        // iso == 2 || iso == -2 (n pbar or p nbar)

        const std::vector<G4double> BFMM185 = {-0.734, 0.841, 0.905, 3.415, -2.316, 0.775};
        //{22.781, -22.602, -0.752, -11.036, 1.548, 0.775}
        //const G4double Eth_PPbar_PPbar_pi0 = 0.775;
        const std::vector<G4double> BFMM188 = { -0.442, 0.501, 0.002, 3.434, -1.201, 0.798};
        //const G4double Eth_PPbar_PNbar_pim = 0.798;
        const std::vector<G4double> BFMM199 = {-2.025, 2.055, -2.355, 6.064, -2.004, 0.798};
        //const G4double Eth_PPbar_NPbar_pip = 0.798;
        const std::vector<G4double> BFMM491 = {24.125, -20.669, -1.534, -19.573, 4.493, 0.787};
        //const G4double Eth_NPbar_PPbar_pim = 0.787; 
        const std::vector<G4double> BFMM495 = {-0.650, -0.140, -0.058, 5.166, -1.705, 0.777};
        //const G4double Eth_NPbar_NPbar_pi0 = 0.777;

        const Particle *antinucleon;
        const Particle *nucleon;
        
        if (p1->isAntiNucleon()) {
            antinucleon = p1;
            nucleon = p2;
        }
        else {
            antinucleon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antinucleon, nucleon); // GeV
        
        if(iso == 2 || iso == -2){ //npbar or pnbar
            sigma = KinematicsUtils::compute_xs(BFMM491, pLab) + KinematicsUtils::compute_xs(BFMM185, pLab) + KinematicsUtils::compute_xs(BFMM188, pLab);
            return sigma;
        }
        else{ // ppbar or nnbar
            sigma = KinematicsUtils::compute_xs(BFMM199, pLab) + KinematicsUtils::compute_xs(BFMM185, pLab) + KinematicsUtils::compute_xs(BFMM188, pLab);
            return sigma;
        }
    }

    G4double CrossSectionsAntiparticles::NNbarToNNbar2pi(Particle const * const p1, Particle const * const p2) {
        //brief ppbar
        // p pbar -> p pbar pi+ pi- (BFMM 167)
        // p pbar -> p nbar pi- pi0 (same as BFMM 490)
        // p pbar -> n pbar pi+ pi0 (same as BFMM 490)
        // p pbar -> n nbar pi+ pi- (BFMM 198)
        //
        //brief npbar
        // n pbar -> p pbar pi- pi0 (BFMM 490)
        // n pbar -> p nbar pi- pi- (BFMM 492)
        // n pbar -> n pbar pi+ pi- (BFMM 494)
        // n pbar -> n nbar pi- pi0 (same as BFMM 490)
        //
        //brief nnbar
        // n nbar -> n nbar pi+ pi- (same as BFMM 167)
        // n nbar -> p nbar pi- pi0 (same as BFMM 490)
        // n nbar -> n pbar pi+ pi0 (same as BFMM 490)
        // n nbar -> p pbar pi+ pi- (same as BFMM 198)
        //
        //brief pnbar
        // p nbar -> p pbar pi+ pi0 (same as BFMM 490)
        // p nbar -> n pbar pi+ pi+ (same as BFMM 492)
        // p nbar -> p nbar pi+ pi- (same as BFMM 494)
        // p nbar -> n nbar pi+ pi0 (same as BFMM 490)
        //
        //
        // BFMM 188,199 are very close in value, 491 is larger

// assert((p1->isAntiNucleon() && p2->isNucleon()) || (p1->isNucleon() && p2->isAntiNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        // iso == 2 || iso == -2 (n pbar or p nbar)

        const std::vector<G4double> BFMM167 = {-6.885, 0.476, 1.206, 13.857, -5.728, 1.220};
        //const G4double Eth_PPbar_PPbar_pip_pim = 1.220;
        const std::vector<G4double> BFMM198 = {1.857, -21.213, -3.448, 0.827, -0.390, 1.231};
        //const G4double Eth_PPbar_NNbar_pip_pim = 1.231;
        const std::vector<G4double> BFMM490 = {-3.594, 0.811, 0.306, 5.108, -1.625, 1.201};
        //const G4double Eth_PNbar_PPbar_pim_pi0 = 1.201;
        const std::vector<G4double> BFMM492 = {-5.443, 7.254, -2.936, 8.441, -2.588, 1.221};
        //const G4double Eth_PNbar_NPbar_pim_pim = 1.221;
        const std::vector<G4double> BFMM494 = {21.688, -38.709, -2.062, -17.783, 3.895, 1.221};
        //const G4double Eth_NPbar_NPbar_pip_pim = 1.221; 

        const Particle *antinucleon;
        const Particle *nucleon;
        
        if (p1->isAntiNucleon()) {
            antinucleon = p1;
            nucleon = p2;
        }
        else {
            antinucleon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antinucleon, nucleon); // GeV
        
        if(iso == 2 || iso == -2){ // pnbar or npbar
            sigma = KinematicsUtils::compute_xs(BFMM490, pLab) + KinematicsUtils::compute_xs(BFMM490, pLab) + KinematicsUtils::compute_xs(BFMM167, pLab) + KinematicsUtils::compute_xs(BFMM198, pLab);
            return sigma;
        }
        else{ // ppbar or nnbar
            sigma = KinematicsUtils::compute_xs(BFMM490, pLab) + KinematicsUtils::compute_xs(BFMM490, pLab) + KinematicsUtils::compute_xs(BFMM492, pLab) + KinematicsUtils::compute_xs(BFMM494, pLab);
            return sigma;
        }
    }

    G4double CrossSectionsAntiparticles::NNbarToNNbar3pi(Particle const * const p1, Particle const * const p2) {
        //brief ppbar
        // p pbar -> p pbar pi+ pi- pi0 (BFMM 161)
        // p pbar -> p nbar 2pi- pi+ (BFMM 169)
        // p pbar -> n pbar 2pi+ pi- (BFMM 201)
        // p pbar -> n nbar pi+ pi- pi0 (BFMM 197)
        //
        //brief npbar
        // n pbar -> p pbar 2pi- pi+ (same as BFMM 169)
        // n pbar -> p nbar 2pi- pi0 (same as BFMM 197)
        // n pbar -> n pbar pi+ pi- pi0 (same as BFMM 161)
        // n pbar -> n nbar 2pi- pi+ (same as BFMM 169)
        //
        //brief nnbar
        // n nbar -> n nbar pi+ pi- pi0 (same as BFMM 161)
        // n nbar -> p nbar 2pi- pi+ (same as BFMM 169)
        // n nbar -> n pbar 2pi+ pi- (same as BFMM 201)
        // n nbar -> p pbar pi+ pi- pi0 (same as BFMM 197)
        //
        //brief pnbar
        // p nbar -> p pbar 2pi+ pi- (same as BFMM 169)
        // p nbar -> n pbar 2pi+ pi0 (same as BFMM 197)
        // p nbar -> p nbar pi+ pi- pi0 (same as BFMM 161)
        // p nbar -> n nbar 2pi+ pi- (same as BFMM 169)
        //      

// assert((p1->isAntiNucleon() && p2->isNucleon()) || (p1->isNucleon() && p2->isAntiNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        // iso == 2 || iso == -2 (n pbar or p nbar)

        const std::vector<G4double> BFMM161 = {-6.434, 1.351, -5.185, 7.754, -1.692, 1.604};
        //const G4double Eth_PPbar_PPbar_pip_pim_pi0 = 1.604;
        const std::vector<G4double> BFMM169 = {3.696, -5.356, -0.053, 1.941, -0.432, 1.624};
        //const G4double Eth_PPbar_PNbar_2pim_pip = 1.624;
        const std::vector<G4double> BFMM201 = {-1.070, -0.636, -0.009, 2.335, -0.499, 1.624};
        //const G4double Eth_PPbar_NPbar_2pip_pim = 1.624;
        const std::vector<G4double> BFMM197 = {1.857, -21.213, -3.448, 0.827, -0.390, 1.616};
        //const G4double Eth_PPbar_NNbar_pip_pim_pi0 = 1.616;

        const Particle *antinucleon;
        const Particle *nucleon;
        
        if (p1->isAntiNucleon()) {
            antinucleon = p1;
            nucleon = p2;
        }
        else {
            antinucleon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antinucleon, nucleon); // GeV
        
        if(iso == 2 || iso == -2){ // pnbar or npbar
            sigma = KinematicsUtils::compute_xs(BFMM169, pLab) + KinematicsUtils::compute_xs(BFMM169, pLab) + KinematicsUtils::compute_xs(BFMM197, pLab) + KinematicsUtils::compute_xs(BFMM161, pLab);
            return sigma;
        }
        else{ // ppbar or nnbar
            sigma = KinematicsUtils::compute_xs(BFMM161, pLab) + KinematicsUtils::compute_xs(BFMM169, pLab) + KinematicsUtils::compute_xs(BFMM197, pLab) + KinematicsUtils::compute_xs(BFMM201, pLab);
            return sigma;
        }
    }

    G4double CrossSectionsAntiparticles::NNbarToAnnihilation(Particle const * const p1, Particle const * const p2) {
        //brief ppbar
        /*
        This part only contains total annihilation xs, the choice of a particular final state
        will be done in the channel file.
        As long as we only have good data for ppbar, we assume that for npbar, pnbar and nnbar the xs
        will be the same, but in order to compensate for the Coulombic effect the ppbar annihilation xs
        is multiplied by the pnbar total xs and divided by the ppbar total xs.
        */

// assert((p1->isAntiNucleon() && p2->isNucleon()) || (p1->isNucleon() && p2->isAntiNucleon()));
        
        G4double sigma=0.;
        const G4int iso=ParticleTable::getIsospin(p1->getType()) + ParticleTable::getIsospin(p2->getType());
        // iso == 2 || iso == -2 (n pbar or p nbar)

        const std::vector<G4double> BFMM6 = {66.098, 0.153, -4.576, -38.319, 6.625}; //ppbar annihilation xs
        const std::vector<G4double> BFMM1 = {119.066, 6.251, -0.006, -60.046, 11.958}; //ppbar total xs
        const std::vector<G4double> BFMM471 = {108.104, 15.708, 0.832, -54.632, -6.958}; //npbar total xs

        const Particle *antinucleon;
        const Particle *nucleon;
        
        if (p1->isAntiNucleon()) {
            antinucleon = p1;
            nucleon = p2;
        }
        else {
            antinucleon = p2;
            nucleon = p1;
        }
        
        const G4double pLab = 0.001*KinematicsUtils::momentumInLab(antinucleon, nucleon); // GeV
        
        if(iso == 2 || iso == -2){ // pnbar or npbar
            sigma = KinematicsUtils::compute_xs(BFMM6, pLab)*KinematicsUtils::compute_xs(BFMM471, pLab)/KinematicsUtils::compute_xs(BFMM1, pLab);
            return sigma;
        }
        else if(p1->getType()==antiProton || p2->getType()==Proton){ // ppbar case
            sigma = KinematicsUtils::compute_xs(BFMM6, pLab);
            return sigma;
        }
        else{ // nnbar case
            sigma = KinematicsUtils::compute_xs(BFMM6, pLab)*KinematicsUtils::compute_xs(BFMM471, pLab)/KinematicsUtils::compute_xs(BFMM1, pLab);
            return sigma;
        }
    }

} // namespace G4INCL

