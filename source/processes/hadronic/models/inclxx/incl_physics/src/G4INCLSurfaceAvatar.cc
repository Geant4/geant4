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

/*
 * G4INCLReflectionAvatar.cc
 *
 *  \date Jun 8, 2009
 * \author Pekka Kaitaniemi
 */

#include "G4INCLSurfaceAvatar.hh"
#include "G4INCLRandom.hh"
#include "G4INCLReflectionChannel.hh"
#include "G4INCLTransmissionChannel.hh"
#include "G4INCLClustering.hh"
#include <sstream>
#include <string>

namespace G4INCL {

  SurfaceAvatar::SurfaceAvatar(G4INCL::Particle *aParticle, G4double time, G4INCL::Nucleus *n)
    :IAvatar(time), theParticle(aParticle), theNucleus(n),
    particlePIn(0.),
    particlePOut(0.),
    particleTOut(0.),
    TMinusV(0.),
    TMinusV2(0.),
    particleMass(0.),
    sinIncidentAngle(0.),
    cosIncidentAngle(0.),
    sinRefractionAngle(0.),
    cosRefractionAngle(0.),
    refractionIndexRatio(0.),
    internalReflection(false)
  {
    setType(SurfaceAvatarType);
  }

  SurfaceAvatar::~SurfaceAvatar()
  {
  }

  G4INCL::IChannel* SurfaceAvatar::getChannel() {
    if(theParticle->isTargetSpectator()) {
      INCL_DEBUG("Particle " << theParticle->getID() << " is a spectator, reflection" << '\n');
      return new ReflectionChannel(theNucleus, theParticle);
    }

    // We forbid transmission of resonances below the Fermi energy. Emitting a
    // delta particle below Tf can lead to negative excitation energies, since
    // CDPP assumes that particles stay in the Fermi sea.
    if(theParticle->isResonance()) {
      const G4double theFermiEnergy = theNucleus->getPotential()->getFermiEnergy(theParticle);
      if(theParticle->getKineticEnergy()<theFermiEnergy) {
        INCL_DEBUG("Particle " << theParticle->getID() << " is a resonance below Tf, reflection" << '\n'
            << "  Tf=" << theFermiEnergy << ", EKin=" << theParticle->getKineticEnergy() << '\n');
        return new ReflectionChannel(theNucleus, theParticle);
      }
    }

    // Don't try to make a cluster if the leading particle is too slow
    const G4double transmissionProbability = getTransmissionProbability(theParticle);
    const G4double TOut = TMinusV;
    const G4double kOut = particlePOut;
    const G4double cosR = cosRefractionAngle;

    INCL_DEBUG("Transmission probability for particle " << theParticle->getID() << " = " << transmissionProbability << '\n');
    /* Don't attempt to construct clusters when a projectile spectator is
     * trying to escape during a nucleus-nucleus collision. The idea behind
     * this is that projectile spectators will later be collected in the
     * projectile remnant, and trying to clusterise them somewhat feels like
     * G4double counting. Moreover, applying the clustering algorithm on escaping
     * projectile spectators makes the code *really* slow if the projectile is
     * large.
     */
    if(theParticle->isNucleonorLambda()
        && (!theParticle->isProjectileSpectator() || !theNucleus->isNucleusNucleusCollision())
        && transmissionProbability>1.E-4) {
      Cluster *candidateCluster = 0;

      candidateCluster = Clustering::getCluster(theNucleus, theParticle);
      if(candidateCluster != 0 &&
          Clustering::clusterCanEscape(theNucleus, candidateCluster)) {

        INCL_DEBUG("Cluster algorithm succeded. Candidate cluster:" << '\n' << candidateCluster->print() << '\n');

        // Check if the cluster can penetrate the Coulomb barrier
        const G4double clusterTransmissionProbability = getTransmissionProbability(candidateCluster);
        const G4double x = Random::shoot();

        INCL_DEBUG("Transmission probability for cluster " << candidateCluster->getID() << " = " << clusterTransmissionProbability << '\n');

        if (x <= clusterTransmissionProbability) {
          theNucleus->getStore()->getBook().incrementEmittedClusters();
          INCL_DEBUG("Cluster " << candidateCluster->getID() << " passes the Coulomb barrier, transmitting." << '\n');
          return new TransmissionChannel(theNucleus, candidateCluster);
        } else {
          INCL_DEBUG("Cluster " << candidateCluster->getID() << " does not pass the Coulomb barrier. Falling back to transmission of the leading particle." << '\n');
          delete candidateCluster;
        }
      } else {
        delete candidateCluster;
      }
    }

    // If we haven't transmitted a cluster (maybe cluster feature was
    // disabled or maybe we just can't produce an acceptable cluster):

    // Always transmit projectile spectators if no cluster was formed and if
    // transmission is energetically allowed
    if(theParticle->isProjectileSpectator() && transmissionProbability>0.) {
      INCL_DEBUG("Particle " << theParticle->getID() << " is a projectile spectator, transmission" << '\n');
      return new TransmissionChannel(theNucleus, theParticle, TOut);
    }

    // Transmit or reflect depending on the transmission probability
    const G4double x = Random::shoot();

    if(x <= transmissionProbability) { // Transmission
      INCL_DEBUG("Particle " << theParticle->getID() << " passes the Coulomb barrier, transmitting." << '\n');
      if(theNucleus->getStore()->getConfig()->getRefraction()) {
        return new TransmissionChannel(theNucleus, theParticle, kOut, cosR);
      } else {
        return new TransmissionChannel(theNucleus, theParticle, TOut);
      }
    } else { // Reflection
      INCL_DEBUG("Particle " << theParticle->getID() << " does not pass the Coulomb barrier, reflection." << '\n');
      return new ReflectionChannel(theNucleus, theParticle);
    }
  }

  void SurfaceAvatar::fillFinalState(FinalState *fs) {
    getChannel()->fillFinalState(fs);
  }

  void SurfaceAvatar::preInteraction() {}

  void SurfaceAvatar::postInteraction(FinalState *fs) {
    ParticleList const &outgoing = fs->getOutgoingParticles();
    if(!outgoing.empty()) { // Transmission
// assert(outgoing.size()==1);
      Particle *out = outgoing.front();
      out->rpCorrelate();
      if(out->isCluster()) {
        Cluster *clusterOut = dynamic_cast<Cluster*>(out);
        ParticleList const &components = clusterOut->getParticles();
        for(ParticleIter i=components.begin(), e=components.end(); i!=e; ++i) {
          if(!(*i)->isTargetSpectator())
            theNucleus->getStore()->getBook().decrementCascading();
        }
        out->setBiasCollisionVector(components.getParticleListBiasVector());
      } else if(!theParticle->isTargetSpectator()) {
// assert(out==theParticle);
        theNucleus->getStore()->getBook().decrementCascading();
      }
    }
  }

  std::string SurfaceAvatar::dump() const {
    std::stringstream ss;
    ss << "(avatar " << theTime << " 'reflection" << '\n'
      << "(list " << '\n'
      << theParticle->dump()
      << "))" << '\n';
    return ss.str();
  }

  G4double SurfaceAvatar::getTransmissionProbability(Particle const * const particle) {

    particleMass = particle->getMass();
    const G4double V = particle->getPotentialEnergy();

    // Correction to the particle kinetic energy if using real masses
    const G4int theA = theNucleus->getA();
    const G4int theZ = theNucleus->getZ();
    const G4double correction = particle->getEmissionQValueCorrection(theA, theZ);
    particleTOut = particle->getKineticEnergy() + correction;

    if (particleTOut <= V) // No transmission if total energy < 0
      return 0.0;

    TMinusV = particleTOut-V;
    TMinusV2 = TMinusV*TMinusV;

    // Momenta in and out
    const G4double particlePIn2  = particle->getMomentum().mag2();
    const G4double particlePOut2 = 2.*particleMass*TMinusV+TMinusV2;
    particlePIn  = std::sqrt(particlePIn2);
    particlePOut = std::sqrt(particlePOut2);
    
    if (0. > V) // Automatic transmission for repulsive potential
      return 1.0;

    // Compute the transmission probability
    G4double theTransmissionProbability;
    if(theNucleus->getStore()->getConfig()->getRefraction()) {
      // Use the formula with refraction
      initializeRefractionVariables(particle);

      if(internalReflection)
        return 0.; // total internal reflection

      // Intermediate variables for calculation
      const G4double x = refractionIndexRatio*cosIncidentAngle;
      const G4double y = (x - cosRefractionAngle) / (x + cosRefractionAngle);

      theTransmissionProbability = 1. - y*y;
    } else {
      // Use the formula without refraction
      // Intermediate variable for calculation
      const G4double y = particlePIn+particlePOut;

      // The transmission probability for a potential step
      theTransmissionProbability = 4.*particlePIn*particlePOut/(y*y);
    }

    // For neutral and negative particles, no Coulomb transmission
    // Also, no Coulomb if the particle takes away all of the nuclear charge
    const G4int particleZ = particle->getZ();
    if (particleZ <= 0 || particleZ >= theZ)
      return theTransmissionProbability;

    // Nominal Coulomb barrier
    const G4double theTransmissionBarrier = theNucleus->getTransmissionBarrier(particle);
    if (TMinusV >= theTransmissionBarrier) // Above the Coulomb barrier
      return theTransmissionProbability;

    // Coulomb-penetration factor
    const G4double px = std::sqrt(TMinusV/theTransmissionBarrier);
    const G4double logCoulombTransmission =
      particleZ*(theZ-particleZ)/137.03*std::sqrt(2.*particleMass/TMinusV/(1.+TMinusV/2./particleMass))
      *(Math::arcCos(px)-px*std::sqrt(1.-px*px));
    INCL_DEBUG("Coulomb barrier, logCoulombTransmission=" << logCoulombTransmission << '\n');
    if (logCoulombTransmission > 35.) // Transmission is forbidden by Coulomb
      return 0.;
    theTransmissionProbability *= std::exp(-2.*logCoulombTransmission);

    return theTransmissionProbability;
  }

  void SurfaceAvatar::initializeRefractionVariables(Particle const * const particle) {
    cosIncidentAngle = particle->getCosRPAngle();
    if(cosIncidentAngle>1.)
      cosIncidentAngle=1.;
    sinIncidentAngle = std::sqrt(1. - cosIncidentAngle*cosIncidentAngle);
    refractionIndexRatio = particlePIn/particlePOut;
    const G4double sinCandidate = refractionIndexRatio*sinIncidentAngle;
    internalReflection = (std::fabs(sinCandidate)>1.);
    if(internalReflection) {
      sinRefractionAngle = 1.;
      cosRefractionAngle = 0.;
    } else {
      sinRefractionAngle = sinCandidate;
      cosRefractionAngle = std::sqrt(1. - sinRefractionAngle*sinRefractionAngle);
    }
    INCL_DEBUG("Refraction parameters initialised as follows:\n"
          << "  cosIncidentAngle=" << cosIncidentAngle << '\n'
          << "  sinIncidentAngle=" << sinIncidentAngle << '\n'
          << "  cosRefractionAngle=" << cosRefractionAngle << '\n'
          << "  sinRefractionAngle=" << sinRefractionAngle << '\n'
          << "  refractionIndexRatio=" << refractionIndexRatio << '\n'
          << "  internalReflection=" << internalReflection << '\n');
  }
}
