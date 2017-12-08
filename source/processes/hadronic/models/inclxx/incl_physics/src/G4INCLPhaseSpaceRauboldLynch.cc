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

#include "G4INCLPhaseSpaceRauboldLynch.hh"
#include "G4INCLRandom.hh"
#include "G4INCLKinematicsUtils.hh"
#include <algorithm>
#include <numeric>
#include <functional>

namespace G4INCL {

  const G4double PhaseSpaceRauboldLynch::wMaxMasslessX[PhaseSpaceRauboldLynch::wMaxNE] = {
    0.01,
    0.017433288222,
    0.0303919538231,
    0.0529831690628,
    0.0923670857187,
    0.161026202756,
    0.280721620394,
    0.489390091848,
    0.853167852417,
    1.48735210729,
    2.5929437974,
    4.52035365636,
    7.88046281567,
    13.7382379588,
    23.9502661999,
    41.7531893656,
    72.7895384398,
    126.896100317,
    221.221629107,
    385.662042116,
    672.33575365,
    1172.10229753,
    2043.35971786,
    3562.24789026,
    6210.16941892,
    10826.3673387,
    18873.9182214,
    32903.4456231,
    57361.5251045,
    100000.0
  };

  const G4double PhaseSpaceRauboldLynch::wMaxMasslessY[PhaseSpaceRauboldLynch::wMaxNE] = {
    -4.69180212056,
    -4.13600582212,
    -3.58020940928,
    -3.02441303018,
    -2.46861663471,
    -1.91282025644,
    -1.35702379603,
    -0.801227342882,
    -0.245430866403,
    0.310365269464,
    0.866161720461,
    1.42195839972,
    1.97775459763,
    2.5335510251,
    3.08934734235,
    3.64514378357,
    4.20094013304,
    4.75673663205,
    5.3125329676,
    5.86832950014,
    6.42412601468,
    6.97992197797,
    7.53571856765,
    8.09151503031,
    8.64731144398,
    9.20310774148,
    9.7589041009,
    10.314700625,
    10.8704971207,
    11.4262935304
  };

  const G4double PhaseSpaceRauboldLynch::wMaxCorrectionX[PhaseSpaceRauboldLynch::wMaxNE] = {
    8.77396538453e-07,
    1.52959067398e-06,
    2.66657950812e-06,
    4.6487249132e-06,
    8.10425612766e-06,
    1.41283832898e-05,
    2.46304178003e-05,
    4.2938917254e-05,
    7.4856652043e-05,
    0.00013049975904,
    0.000227503991225,
    0.000396614265067,
    0.000691429079588,
    0.00120538824295,
    0.00210138806588,
    0.00366341038188,
    0.00638652890627,
    0.0111338199161,
    0.0194099091609,
    0.0338378540766,
    0.0589905062931,
    0.102839849857,
    0.179283674326,
    0.312550396803,
    0.544878115136,
    0.949901722703,
    1.65599105145,
    2.88693692929,
    5.03288035671,
    8.77396538453
  };
  const G4double PhaseSpaceRauboldLynch::wMaxCorrectionY[PhaseSpaceRauboldLynch::wMaxNE] = {
    7.28682764024,
    7.0089298602,
    6.7310326046,
    6.45313436085,
    6.17523713298,
    5.89734080335,
    5.61944580818,
    5.34155326611,
    5.06366530628,
    4.78578514804,
    4.50791742491,
    4.23007259554,
    3.97566018117,
    3.72116551935,
    3.44355099732,
    3.1746329047,
    2.89761210229,
    2.62143335535,
    2.34649707964,
    2.07366126445,
    1.8043078447,
    1.54056992953,
    1.27913996411,
    0.981358535322,
    0.73694629682,
    0.587421056631,
    0.417178268911,
    0.280881750176,
    0.187480311508,
    0.116321254846
  };

  const G4double PhaseSpaceRauboldLynch::wMaxInterpolationMargin = std::log(1.5);

  PhaseSpaceRauboldLynch::PhaseSpaceRauboldLynch() :
    nParticles(0),
    sqrtS(0.),
    availableEnergy(0.),
    maxGeneratedWeight(0.)
  {
    std::vector<G4double> wMaxMasslessXV(wMaxMasslessX, wMaxMasslessX + wMaxNE);
    std::vector<G4double> wMaxMasslessYV(wMaxMasslessY, wMaxMasslessY + wMaxNE);
    wMaxMassless = new InterpolationTable(wMaxMasslessXV, wMaxMasslessYV);
    std::vector<G4double> wMaxCorrectionXV(wMaxCorrectionX, wMaxCorrectionX + wMaxNE);
    std::vector<G4double> wMaxCorrectionYV(wMaxCorrectionY, wMaxCorrectionY + wMaxNE);
    wMaxCorrection = new InterpolationTable(wMaxCorrectionXV, wMaxCorrectionYV);

    // initialize the precalculated coefficients
    prelog[0] = 0.;
    for(size_t i=1; i<wMaxNP; ++i) {
      prelog[i] = -std::log(G4double(i));
    }
  }

  PhaseSpaceRauboldLynch::~PhaseSpaceRauboldLynch() {
    delete wMaxMassless;
    delete wMaxCorrection;
  }

  void PhaseSpaceRauboldLynch::generate(const G4double sqs, ParticleList &particles) {
    maxGeneratedWeight = 0.;

    sqrtS = sqs;

    // initialize the structures containing the particle masses
    initialize(particles);

    // initialize the maximum weight
    const G4double weightMax = computeMaximumWeightParam();

    const G4int maxIter = 500;
    G4int iter = 0;
    G4double weight, r;
    do {
      weight = computeWeight();
      maxGeneratedWeight = std::max(weight, maxGeneratedWeight);
// assert(weight<=weightMax);
      r = Random::shoot();
    } while(++iter<maxIter && r*weightMax>weight); /* Loop checking, 10.07.2015, D.Mancusi */

#ifndef INCLXX_IN_GEANT4_MODE
    if(iter>=maxIter) {
      INCL_WARN("Number of tries exceeded in PhaseSpaceRauboldLynch::generate()\n"
                << "nParticles=" << nParticles << ", sqrtS=" << sqrtS << '\n');
    }
#endif

    generateEvent(particles);
  }

  void PhaseSpaceRauboldLynch::initialize(ParticleList &particles) {
    nParticles = particles.size();
// assert(nParticles>2);

    // masses and sum of masses
    masses.resize(nParticles);
    sumMasses.resize(nParticles);
    std::transform(particles.begin(), particles.end(), masses.begin(), std::mem_fun(&Particle::getMass));
    std::partial_sum(masses.begin(), masses.end(), sumMasses.begin());

    // sanity check
    availableEnergy  = sqrtS-sumMasses[nParticles-1];
// assert(availableEnergy>-1.e-5);
    if(availableEnergy<0.)
      availableEnergy = 0.;

    rnd.resize(nParticles);
    invariantMasses.resize(nParticles);
    momentaCM.resize(nParticles-1);
  }

  G4double PhaseSpaceRauboldLynch::computeMaximumWeightNaive() {
    // compute the maximum weight
    G4double eMMax = sqrtS + masses[0];
    G4double eMMin = 0.;
    G4double wMax = 1.;
    for(size_t i=1; i<nParticles; i++) {
      eMMin += masses[i-1];
      eMMax += masses[i];
      wMax *= KinematicsUtils::momentumInCM(eMMax, eMMin, masses[i]);
    }
    return wMax;
  }

  G4double PhaseSpaceRauboldLynch::computeMaximumWeightParam() {
#ifndef INCLXX_IN_GEANT4_MODE
    if(nParticles>=wMaxNP) {
      INCL_WARN("The requested number of particles (" << nParticles << ") requires extrapolation the tables in PhaseSpaceRauboldLynch. YMMV." << '\n');
    }

    if(availableEnergy>wMaxMasslessX[wMaxNE-1] || availableEnergy<wMaxMasslessX[0] ) {
      INCL_WARN("The requested available energy (" << availableEnergy << " MeV) requires extrapolation the tables in PhaseSpaceRauboldLynch. YMMV." << '\n');
    }
#endif

    // compute the maximum weight
    const G4double logMassless = ((*wMaxMassless)(availableEnergy)+prelog[nParticles])*(nParticles-1);
    const G4double reducedSqrtS = availableEnergy/sumMasses[nParticles-1];
    const G4double correction = (*wMaxCorrection)(reducedSqrtS);
    const G4double wMax = std::exp(correction*G4double(nParticles-1) + logMassless + wMaxInterpolationMargin);
    if(wMax>0.)
      return wMax;
    else
      return computeMaximumWeightNaive();
  }

  G4double PhaseSpaceRauboldLynch::computeWeight() {
    // Generate nParticles-2 sorted random numbers, bracketed by 0 and 1
    rnd[0] = 0.;
    for(size_t i=1; i<nParticles-1; ++i)
      rnd[i] = Random::shoot();
    rnd[nParticles-1] = 1.;
    std::sort(rnd.begin()+1, rnd.begin()+nParticles-1);
    
    // invariant masses
    for(size_t i=0; i<nParticles; ++i)
      invariantMasses[i] = rnd[i]*availableEnergy + sumMasses[i];

    // compute the CM momenta and the weight for the current event
    G4double weight=KinematicsUtils::momentumInCM(invariantMasses[1], invariantMasses[0], masses[1]);
    momentaCM[0] = weight;
    for(size_t i=1; i<nParticles-1; ++i) {
      G4double momentumCM;
// assert(invariantMasses[i+1]-invariantMasses[i]-masses[i+1]> - 1.E-8);
      if(invariantMasses[i+1]-invariantMasses[i]-masses[i+1] < 0.) momentumCM = 0.;
      else momentumCM = KinematicsUtils::momentumInCM(invariantMasses[i+1], invariantMasses[i], masses[i+1]);
      momentaCM[i] = momentumCM;
      weight *= momentumCM;
    }

    return weight;
  }

  void PhaseSpaceRauboldLynch::generateEvent(ParticleList &particles) {
    Particle *p = particles[0];
    ThreeVector momentum = Random::normVector(momentaCM[0]);
    p->setMomentum(momentum);
    p->adjustEnergyFromMomentum();

    ThreeVector boostV;

    for(size_t i=1; i<nParticles; ++i) {
      p = particles[i];
      p->setMomentum(-momentum);
      p->adjustEnergyFromMomentum();

      if(i==nParticles-1)
        break;

      momentum = Random::normVector(momentaCM[i]);

      const G4double iM = invariantMasses[i];
      const G4double recoilE = std::sqrt(momentum.mag2() + iM*iM);
      boostV = -momentum/recoilE;
      for(size_t j=0; j<=i; ++j)
        particles[j]->boost(boostV);
    }
  }

  G4double PhaseSpaceRauboldLynch::getMaxGeneratedWeight() const {
    return maxGeneratedWeight;
  }
}
