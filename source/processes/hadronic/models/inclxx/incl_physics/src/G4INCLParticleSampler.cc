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
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLParticleSampler.cc
 * \brief Class for sampling particles in a nucleus
 *
 * \date 18 July 2012
 * \author Davide Mancusi
 */

#include "G4INCLParticleSampler.hh"

namespace G4INCL {

  ParticleSampler::ParticleSampler(const G4int A, const G4int Z, InverseInterpolationTable const * const rCDFTable, InverseInterpolationTable const * const pCDFTable) :
    sampleOneParticle(&ParticleSampler::sampleOneParticleWithoutRPCorrelation),
    theA(A),
    theZ(Z),
    theRCDFTable(rCDFTable),
    thePCDFTable(pCDFTable),
    theDensity(NULL),
    thePotential(NULL)
  {}

  ParticleSampler::~ParticleSampler() {
  }

  void ParticleSampler::setDensity(NuclearDensity const * const d) {
    theDensity = d;
    updateSampleOneParticleMethod();
  }

  void ParticleSampler::setPotential(NuclearPotential::INuclearPotential const * const p) {
    thePotential = p;
    updateSampleOneParticleMethod();
  }

  void ParticleSampler::updateSampleOneParticleMethod() {
    if(theDensity && thePotential)
      sampleOneParticle = &ParticleSampler::sampleOneParticleWithRPCorrelation;
    else
      sampleOneParticle = &ParticleSampler::sampleOneParticleWithoutRPCorrelation;
  }

  ParticleList ParticleSampler::sampleParticles(const ThreeVector &position) const {
    ParticleList theList;
    if(theA > 2) {
      ParticleType type = Proton;
      for(G4int i = 1; i <= theA; ++i) {
        if(i == (theZ + 1)) // Nucleons [Z+1..A] are neutrons
          type = Neutron;
        Particle *p = (this->*(this->sampleOneParticle))(type);
        p->setPosition(position + p->getPosition());
        theList.push_back(p);
      }
    } else {
      // For deuterons, only sample the proton position and momentum. The
      // neutron position and momenta are determined by the conditions of
      // vanishing CM position and total momentum.
// assert(theZ==1);
      Particle *aProton = (this->*(this->sampleOneParticle))(Proton);
      Particle *aNeutron = new Particle(Neutron, -aProton->getMomentum(), position - aProton->getPosition());
      aProton->setPosition(position + aProton->getPosition());
      theList.push_back(aProton);
      theList.push_back(aNeutron);
    }

    return theList;
  }

  Particle *ParticleSampler::sampleOneParticleWithRPCorrelation(const ParticleType t) const {
// assert(theDensity && thePotential);
    const G4double theFermiMomentum = thePotential->getFermiMomentum(t);
    const ThreeVector momentumVector = Random::sphereVector(theFermiMomentum);
    const G4double momentumRatio = momentumVector.mag()/theFermiMomentum;
    const ThreeVector positionVector = Random::sphereVector(theDensity->getMaxRFromP(momentumRatio));
    return new Particle(t, momentumVector, positionVector);
  }

  Particle *ParticleSampler::sampleOneParticleWithoutRPCorrelation(const ParticleType t) const {
    const G4double position = (*theRCDFTable)(Random::shoot());
    const G4double momentum = (*thePCDFTable)(Random::shoot());
    ThreeVector positionVector = Random::normVector(position);
    ThreeVector momentumVector = Random::normVector(momentum);
    return new Particle(t, momentumVector, positionVector);
  }

}

