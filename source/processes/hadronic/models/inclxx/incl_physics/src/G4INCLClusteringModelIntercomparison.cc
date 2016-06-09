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
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLClusteringModelIntercomparison.hh"
#include "G4INCLCluster.hh"
#include "G4INCLRandom.hh"

namespace G4INCL {
  const G4double ClusteringModelIntercomparison::limitCosEscapeAngle = 0.7;

  static G4bool participantsFirstPredicate(Particle *lhs, Particle * /*rhs*/) {
    return lhs->isParticipant();
  }

  Cluster* ClusteringModelIntercomparison::getCluster(Nucleus *nucleus, Particle *particle) {
    theNucleus = nucleus;
    theLeadingParticle = particle;

    sqtot = 50000.0;
    selectedA = 0;
    selectedZ = 0;

    const G4double transp = 1.0;
    const G4double rmaxws = theNucleus->getDensity()->getMaximumRadius();

    const G4double Rprime = theNucleus->getDensity()->getCentralRadius() + transp;
    const G4double pk = theLeadingParticle->getMomentum().mag();

    const G4double cospr = theLeadingParticle->getPosition().dot(theLeadingParticle->getMomentum())/(theNucleus->getDensity()->getMaximumRadius() * pk);
    const G4double arg = rmaxws*rmaxws - Rprime*Rprime;
    G4double translat = 0.0;

    if(arg > 0.0) {
      const G4double cosmin = std::sqrt(arg)/rmaxws;
      if(cospr <= cosmin) {
	translat = rmaxws * cospr;
      } else {
	translat = rmaxws * (cospr - std::sqrt(cospr*cospr - cosmin*cosmin));
      }
    } else {
      translat = rmaxws * cospr - std::sqrt(Rprime*Rprime - rmaxws*rmaxws*(1.0 - cospr*cospr));
    }

    const ThreeVector oldLeadingParticlePosition = theLeadingParticle->getPosition();
    const ThreeVector leadingParticlePosition = oldLeadingParticlePosition - theLeadingParticle->getMomentum() * (translat/pk);
    const ThreeVector leadingParticleMomentum = theLeadingParticle->getMomentum();
    theLeadingParticle->setPosition(leadingParticlePosition);

    // Select the subset of nucleons that will be considered in the
    // cluster production:
    participantEnergyPool = 0.;
    const ParticleList particles = theNucleus->getStore()->getParticles();
    for(ParticleIter i = particles.begin(); i != particles.end(); ++i) {
      if (!(*i)->isNucleon()) continue; // Only nucleons are allowed in clusters
      if ((*i)->getID() == theLeadingParticle->getID()) continue; // Don't count the leading particle

      G4double space = ((*i)->getPosition() - leadingParticlePosition).mag2();
      G4double momentum = ((*i)->getMomentum() - leadingParticleMomentum).mag2();
      G4double size = space*momentum*ParticleTable::clusterPosFact2[IClusteringModel::maxClusterAlgorithmMass];
      if(size < ParticleTable::clusterPhaseSpaceCut[IClusteringModel::maxClusterAlgorithmMass]) {
	consideredPartners.push_back((*i));
        if((*i)->isParticipant())
          participantEnergyPool += (*i)->getEnergy() - (*i)->getPotentialEnergy() - 931.3;
      }
    }
    // Sort the list of considered partners so that we give priority
    // to participants. As soon as we encounter the first spectator in
    // the list we know that all the remaining nucleons will be
    // spectators too.
    consideredPartners.sort(participantsFirstPredicate);

    runningConfiguration.push_back(theLeadingParticle);
    runningPositions[1] = theLeadingParticle->getPosition();
    runningMomenta[1] = theLeadingParticle->getMomentum();
    runningEnergies[1] = theLeadingParticle->getEnergy();
    runningPotentials[1] = theLeadingParticle->getPotentialEnergy();

    // Start the cluster search!
    findClusterStartingFrom(1, theLeadingParticle->getZ());

    Cluster *chosenCluster = 0;
    if(selectedA!=0) { // A cluster was found!
      chosenCluster =  new Cluster(candidateConfiguration);
    }

    // Restore the original position of the leading particle
    theLeadingParticle->setPosition(oldLeadingParticlePosition);

    cleanUp();
    zeroOut();
    return chosenCluster;
  }

  G4double ClusteringModelIntercomparison::getPhaseSpace(G4int oldA, Particle *p) {
    const G4double psSpace = (p->getPosition() - runningPositions[oldA]).mag2();
    const G4double psMomentum = (p->getMomentum()*oldA - runningMomenta[oldA]).mag2();
    return psSpace * psMomentum * ParticleTable::clusterPosFact2[oldA + 1];
  }

  void ClusteringModelIntercomparison::findClusterStartingFrom(const G4int oldA, const G4int oldZ) {
    const G4int newA = oldA + 1;
    G4int newZ = 0;
    G4int newN = 0;

    for(ParticleIter i = consideredPartners.begin(); i != consideredPartners.end();
	++i) {
      // Only accept particles that are not already part of the cluster
      if((*i)->isInList(runningConfiguration)) continue;

      newZ = oldZ + (*i)->getZ();
      newN = newA - newZ;

      // Skip this nucleon if we already have too many protons or neutrons
      if(newZ > clusterZMaxAll || newN > clusterNMaxAll)
        continue;

      // Compute the phase space factor for a new cluster which
      // consists of the previous running cluster and the new
      // candidate nucleon:
      const G4double phaseSpace = getPhaseSpace(oldA, (*i));
      if(phaseSpace > ParticleTable::clusterPhaseSpaceCut[newA]) continue;

      // eclst:
      runningEnergies[newA] = runningEnergies[oldA] + (*i)->getEnergy();
      // vcl:
      runningPotentials[newA] = runningPotentials[oldA] + (*i)->getPotentialEnergy();

      // Update the available participant kinetic energy
      G4double oldParticipantEnergyPool = participantEnergyPool;
      if((*i)->isParticipant())
        participantEnergyPool -= (*i)->getEnergy() - (*i)->getPotentialEnergy() - 931.3;

      // Check an approximate Coulomb barrier
      const G4double halfB = 0.72 * newZ * theNucleus->getZ()/(theNucleus->getDensity()->getCentralRadius()+1.7);
      const G4double tout = runningEnergies[newA] - runningPotentials[newA] - 931.3*newA;
      if(tout<=halfB && tout+participantEnergyPool<=halfB) {
        participantEnergyPool = oldParticipantEnergyPool;
        continue;
      }

      // Accept the nucleon in the cluster
      runningConfiguration.push_back((*i));
      runningPositions[newA] = (runningPositions[oldA] * oldA + (*i)->getPosition())*ParticleTable::clusterPosFact[newA];
      runningMomenta[newA] = runningMomenta[oldA] + (*i)->getMomentum();

      // Keep track of the best physical cluster
      if(newZ >= ParticleTable::clusterZMin[newA] && newZ <= ParticleTable::clusterZMax[newA]) {
        // Note: sqc is real kinetic energy, not the square of the kinetic energy!
        G4double sqc = KinematicsUtils::invariantMass(runningEnergies[newA], runningMomenta[newA]);
        G4double sqct = (sqc - newZ * 938.27
            - (newA - newZ) * 939.57
            - ParticleTable::binding[newZ][newA])
          *ParticleTable::clusterPosFact[newA];

        if(sqct < sqtot) {
          sqtot = sqct;
          selectedA = newA;
          selectedZ = newZ;
          delete candidateConfiguration;
          candidateConfiguration = new ParticleList(runningConfiguration);
        }
      }

      if(newA < IClusteringModel::maxClusterAlgorithmMass && newA+1 < theNucleus->getA()) {
	findClusterStartingFrom(newA, newZ);
      }

      runningConfiguration.pop_back();
      participantEnergyPool = oldParticipantEnergyPool;
    }
  }

  G4bool ClusteringModelIntercomparison::clusterCanEscape(Cluster const * const c) {
    // Check the escape angle of the cluster
    const ThreeVector &pos = c->getPosition();
    const ThreeVector &mom = c->getMomentum();
    const G4double cosEscapeAngle = pos.dot(mom) / std::sqrt(pos.mag2()*mom.mag2());
    if(cosEscapeAngle < limitCosEscapeAngle)
      return false;

    // Check if the cluster can penetrate the Coulomb barrier
    const G4double transmissionProbability = theNucleus->getTransmissionProbability(c);
    const G4double x = Random::shoot();

    return (x <= transmissionProbability);
  }

}
