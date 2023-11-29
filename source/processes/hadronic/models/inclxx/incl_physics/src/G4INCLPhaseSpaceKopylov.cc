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

#include "G4INCLPhaseSpaceKopylov.hh"
#include "G4INCLRandom.hh"
#include "G4INCLKinematicsUtils.hh"
#include <algorithm>
#include <numeric>
#include <functional>

namespace G4INCL {

  G4double PhaseSpaceKopylov::betaKopylov(G4int K) const {
    G4int N = 3*K - 5;
    G4double xN = G4double(N);
    G4double Fmax = std::sqrt(std::pow(xN/(xN+1.),N)/(xN+1.));

    G4double F, chi;
    unsigned long loopCounter = 0;
    const unsigned long maxLoopCounter = 10000000;
    do {
      chi = Random::shoot();
      F = std::sqrt(std::pow(chi,N)*(1.-chi));
      ++loopCounter;
    } while (loopCounter<maxLoopCounter && Fmax*Random::shoot() > F); /* Loop checking, 10.07.2015, D.Mancusi */
    return chi;
  }

  void PhaseSpaceKopylov::generate(const G4double sqrtS, ParticleList &particles) {

    boostV.setX(0.0);
    boostV.setY(0.0);
    boostV.setZ(0.0);

    const std::size_t N = particles.size();
    masses.resize(N);
    sumMasses.resize(N);
    std::transform(particles.begin(), particles.end(), masses.begin(), std::mem_fn(&Particle::getMass));
    std::partial_sum(masses.begin(), masses.end(), sumMasses.begin());

    G4double PFragMagCM = 0.0;
    G4double T = sqrtS-sumMasses.back();
// assert(T>-1.e-5);
    if(T<0.)
      T=0.;

    // The first particle in the list will pick up all the recoil
    Particle *restParticle = particles.front();
    restParticle->setMass(sqrtS);
    restParticle->adjustEnergyFromMomentum();

    G4int k=G4int(N-1);
    for (auto p=particles.rbegin(); k>0; ++p, --k) {
      const G4double mu = sumMasses[k-1];
      T *= (k>1) ? betaKopylov(k) : 0.;

      const G4double restMass = mu + T;

      PFragMagCM = KinematicsUtils::momentumInCM(restParticle->getMass(), masses[k], restMass);
      PFragCM = Random::normVector(PFragMagCM);
      (*p)->setMomentum(PFragCM);
      (*p)->adjustEnergyFromMomentum();
      restParticle->setMass(restMass);
      restParticle->setMomentum(-PFragCM);
      restParticle->adjustEnergyFromMomentum();

      (*p)->boost(boostV);
      restParticle->boost(boostV);

      boostV = -restParticle->boostVector();
    }
    restParticle->setMass(masses[0]);
    restParticle->adjustEnergyFromMomentum();
  }

}
