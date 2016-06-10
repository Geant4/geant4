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

/** \file G4INCLParticleSampler.cc
 * \brief Class for sampling particles in a nucleus
 *
 * \date 18 July 2012
 * \author Davide Mancusi
 */

#include "G4INCLParticleSampler.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLNuclearDensityFactory.hh"

namespace G4INCL {

  ParticleSampler::ParticleSampler(const G4int A, const G4int Z) :
    sampleOneProton(&ParticleSampler::sampleOneParticleWithoutRPCorrelation),
    sampleOneNeutron(&ParticleSampler::sampleOneParticleWithoutRPCorrelation),
    theA(A),
    theZ(Z),
    theDensity(NULL),
    thePotential(NULL)
  {
    std::fill(theRCDFTable, theRCDFTable + UnknownParticle, static_cast<InterpolationTable *>(NULL));
    std::fill(thePCDFTable, thePCDFTable + UnknownParticle, static_cast<InterpolationTable *>(NULL));
    std::fill(rpCorrelationCoefficient, rpCorrelationCoefficient + UnknownParticle, 1.);
    rpCorrelationCoefficient[Proton] = ParticleTable::getRPCorrelationCoefficient(Proton);
    rpCorrelationCoefficient[Neutron] = ParticleTable::getRPCorrelationCoefficient(Neutron);
  }

  ParticleSampler::~ParticleSampler() {
  }

  void ParticleSampler::setDensity(NuclearDensity const * const d) {
    theDensity = d;
    updateSampleOneParticleMethods();
  }

  void ParticleSampler::setPotential(NuclearPotential::INuclearPotential const * const p) {
    thePotential = p;
    updateSampleOneParticleMethods();
  }

  void ParticleSampler::updateSampleOneParticleMethods() {
    if(theDensity && thePotential) {
      if(rpCorrelationCoefficient[Proton]>0.99999) {
        sampleOneProton = &ParticleSampler::sampleOneParticleWithRPCorrelation;
      } else {
        sampleOneProton = &ParticleSampler::sampleOneParticleWithFuzzyRPCorrelation;
      }
      if(rpCorrelationCoefficient[Neutron]>0.99999) {
        sampleOneNeutron = &ParticleSampler::sampleOneParticleWithRPCorrelation;
      } else {
        sampleOneNeutron = &ParticleSampler::sampleOneParticleWithFuzzyRPCorrelation;
      }
    } else {
      sampleOneProton = &ParticleSampler::sampleOneParticleWithoutRPCorrelation;
      sampleOneNeutron = &ParticleSampler::sampleOneParticleWithoutRPCorrelation;
    }
  }

  ParticleList ParticleSampler::sampleParticles(const ThreeVector &position) {
    ParticleList aList;
    sampleParticlesIntoList(position, aList);
    return aList;
  }

  void ParticleSampler::sampleParticlesIntoList(const ThreeVector &position, ParticleList &theList) {

    if(sampleOneProton == &ParticleSampler::sampleOneParticleWithoutRPCorrelation) {
      // sampling without correlation, we need to initialize the CDF tables
      theRCDFTable[Proton] = NuclearDensityFactory::createRCDFTable(Proton, theA, theZ);
      thePCDFTable[Proton] = NuclearDensityFactory::createPCDFTable(Proton, theA, theZ);
      theRCDFTable[Neutron] = NuclearDensityFactory::createRCDFTable(Neutron, theA, theZ);
      thePCDFTable[Neutron] = NuclearDensityFactory::createPCDFTable(Neutron, theA, theZ);
    }

    theList.resize(theA);
    if(theA > 2) {
      ParticleType type = Proton;
      ParticleSamplerMethod sampleOneParticle = sampleOneProton;
      for(G4int i = 0; i < theA; ++i) {
        if(i == theZ) { // Nucleons [Z..A-1] are neutrons
          type = Neutron;
          sampleOneParticle = sampleOneNeutron;
        }
        Particle *p = (this->*sampleOneParticle)(type);
        p->setPosition(position + p->getPosition());
        theList[i] = p;
      }
    } else {
      // For deuterons, only sample the proton position and momentum. The
      // neutron position and momenta are determined by the conditions of
      // vanishing CM position and total momentum.
// assert(theZ==1);
      Particle *aProton = (this->*(this->sampleOneProton))(Proton);
      Particle *aNeutron = new Particle(Neutron, -aProton->getMomentum(), position - aProton->getPosition());
      aProton->setPosition(position + aProton->getPosition());
      theList[0] = aProton;
      theList[1] = aNeutron;
    }
  }

  Particle *ParticleSampler::sampleOneParticleWithRPCorrelation(const ParticleType t) const {
// assert(theDensity && thePotential);
    const G4double theFermiMomentum = thePotential->getFermiMomentum(t);
    const ThreeVector momentumVector = Random::sphereVector(theFermiMomentum);
    const G4double momentumAbs = momentumVector.mag();
    const G4double momentumRatio = momentumAbs/theFermiMomentum;
    const G4double reflectionRadius = theDensity->getMaxRFromP(t, momentumRatio);
    const ThreeVector positionVector = Random::sphereVector(reflectionRadius);
    Particle *aParticle = new Particle(t, momentumVector, positionVector);
    aParticle->setUncorrelatedMomentum(momentumAbs);
    return aParticle;
  }

  Particle *ParticleSampler::sampleOneParticleWithoutRPCorrelation(const ParticleType t) const {
    const G4double position = (*(theRCDFTable[t]))(Random::shoot());
    const G4double momentum = (*(thePCDFTable[t]))(Random::shoot());
    ThreeVector positionVector = Random::normVector(position);
    ThreeVector momentumVector = Random::normVector(momentum);
    return new Particle(t, momentumVector, positionVector);
  }

  Particle *ParticleSampler::sampleOneParticleWithFuzzyRPCorrelation(const ParticleType t) const {
// assert(theDensity && thePotential);
    std::pair<G4double,G4double> ranNumbers = Random::correlatedUniform(rpCorrelationCoefficient[t]);
    const G4double x = Math::pow13(ranNumbers.first);
    const G4double y = Math::pow13(ranNumbers.second);
    const G4double theFermiMomentum = thePotential->getFermiMomentum(t);
    const ThreeVector momentumVector = Random::normVector(y*theFermiMomentum);
    const G4double reflectionRadius = theDensity->getMaxRFromP(t, x);
    const ThreeVector positionVector = Random::sphereVector(reflectionRadius);
    Particle *aParticle = new Particle(t, momentumVector, positionVector);
    aParticle->setUncorrelatedMomentum(x*theFermiMomentum);
    return aParticle;
  }

}

