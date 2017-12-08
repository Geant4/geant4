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

#include "G4INCLClusteringModelIntercomparison.hh"
#include "G4INCLCluster.hh"
#include "G4INCLRandom.hh"
#include "G4INCLHashing.hh"
#include <algorithm>

namespace G4INCL {

  const G4int ClusteringModelIntercomparison::clusterZMin[ParticleTable::maxClusterMass+1] = {0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3};
  const G4int ClusteringModelIntercomparison::clusterZMax[ParticleTable::maxClusterMass+1] = {0, 0, 1, 2, 3, 3, 5, 5, 6, 6, 7, 7, 8};

  const G4double ClusteringModelIntercomparison::clusterPosFact[ParticleTable::maxClusterMass+1] = {0.0, 1.0, 0.5,
    0.33333, 0.25,
    0.2, 0.16667,
    0.14286, 0.125,
    0.11111, 0.1,
    0.09091, 0.083333};

  const G4double ClusteringModelIntercomparison::clusterPosFact2[ParticleTable::maxClusterMass+1] = {0.0, 1.0, 0.25,
    0.11111, 0.0625,
    0.04, 0.0277778,
    0.020408, 0.015625,
    0.012346, 0.01,
    0.0082645, 0.0069444};

  const G4double ClusteringModelIntercomparison::clusterPhaseSpaceCut[ParticleTable::maxClusterMass+1] = {0.0, 70000.0, 180000.0,
    90000.0, 90000.0,
    128941.0 ,145607.0,
    161365.0, 176389.0,
    190798.0, 204681.0,
    218109.0, 231135.0};

  const G4double ClusteringModelIntercomparison::limitCosEscapeAngle = 0.7;

#ifdef INCL_DO_NOT_COMPILE
  namespace {
    G4bool cascadingFirstPredicate(ConsideredPartner const &aPartner) {
      return !aPartner.isTargetSpectator;
    }
  }
#endif

  Cluster* ClusteringModelIntercomparison::getCluster(Nucleus *nucleus, Particle *particle) {
    // Set the maximum clustering mass dynamically, based on the current nucleus
    const G4int maxClusterAlgorithmMass = nucleus->getStore()->getConfig()->getClusterMaxMass();
    runningMaxClusterAlgorithmMass = std::min(maxClusterAlgorithmMass, nucleus->getA()/2);

    // Nucleus too small?
    if(runningMaxClusterAlgorithmMass<=1)
      return NULL;

    theNucleus = nucleus;
    Particle *theLeadingParticle = particle;

    // Initialise sqtot to a large number
    sqtot = 50000.0;
    selectedA = 0;
    selectedZ = 0;

    // The distance parameter, known as h in publications.
    // Default value is 1 fm.
    const G4double transp = 1.0;

    const G4double rmaxws = theNucleus->getUniverseRadius();

    // Radius of the sphere where the leading particle is positioned.
    const G4double Rprime = theNucleus->getDensity()->getProtonNuclearRadius() + transp;

    // Bring the leading particle back to the coalescence sphere
    const G4double pk = theLeadingParticle->getMomentum().mag();
    const G4double cospr = theLeadingParticle->getPosition().dot(theLeadingParticle->getMomentum())/(theNucleus->getUniverseRadius() * pk);
    const G4double arg = rmaxws*rmaxws - Rprime*Rprime;
    G4double translat;

    if(arg > 0.0) {
      // coalescence sphere smaller than Rmax
      const G4double cosmin = std::sqrt(arg)/rmaxws;
      if(cospr <= cosmin) {
        // there is an intersection with the coalescence sphere
        translat = rmaxws * cospr;
      } else {
        // no intersection with the coalescence sphere
        translat = rmaxws * (cospr - std::sqrt(cospr*cospr - cosmin*cosmin));
      }
    } else {
      // coalescence sphere larger than Rmax
      translat = rmaxws * cospr - std::sqrt(Rprime*Rprime - rmaxws*rmaxws*(1.0 - cospr*cospr));
    }

    const ThreeVector oldLeadingParticlePosition = theLeadingParticle->getPosition();
    const ThreeVector leadingParticlePosition = oldLeadingParticlePosition - theLeadingParticle->getMomentum() * (translat/pk);
    const ThreeVector &leadingParticleMomentum = theLeadingParticle->getMomentum();
    theLeadingParticle->setPosition(leadingParticlePosition);

    // Initialise the array of considered nucleons
    const G4int theNucleusA = theNucleus->getA();
    if(nConsideredMax < theNucleusA) {
      delete [] consideredPartners;
      delete [] isInRunningConfiguration;
      nConsideredMax = 2*theNucleusA;
      consideredPartners = new ConsideredPartner[nConsideredMax];
      isInRunningConfiguration = new G4bool [nConsideredMax];
      std::fill(isInRunningConfiguration,
                isInRunningConfiguration + nConsideredMax,
                false);
    }

    // Select the subset of nucleons that will be considered in the
    // cluster production:
    cascadingEnergyPool = 0.;
    nConsidered = 0;
    ParticleList const &particles = theNucleus->getStore()->getParticles();
    for(ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
      if (!(*i)->isNucleonorLambda()) continue; // Only nucleons and lambdas are allowed in clusters
      if ((*i)->getID() == theLeadingParticle->getID()) continue; // Don't count the leading particle

      G4double space = ((*i)->getPosition() - leadingParticlePosition).mag2();
      G4double momentum = ((*i)->getMomentum() - leadingParticleMomentum).mag2();
      G4double size = space*momentum*clusterPosFact2[runningMaxClusterAlgorithmMass];
      // Nucleons are accepted only if they are "close enough" in phase space
      // to the leading nucleon. The selected phase-space parameter corresponds
      // to the running maximum cluster mass.
      if(size < clusterPhaseSpaceCut[runningMaxClusterAlgorithmMass]) {
	consideredPartners[nConsidered] = *i;
        // Keep trace of how much energy is carried by cascading nucleons. This
        // is used to stop the clustering algorithm as soon as possible.
        if(!(*i)->isTargetSpectator())
          cascadingEnergyPool += consideredPartners[nConsidered].energy - consideredPartners[nConsidered].potentialEnergy - 931.3;
        nConsidered++;
        // Make sure we don't exceed the array size
// assert(nConsidered<=nConsideredMax);
      }
    }
    // Sort the list of considered partners so that we give priority
    // to participants. As soon as we encounter the first spectator in
    // the list we know that all the remaining nucleons will be
    // spectators too.
//    std::partition(consideredPartners, consideredPartners+nConsidered, cascadingFirstPredicate);

#ifndef INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_None
    // Clear the sets of checked configurations
    // We stop caching two masses short of the max mass -- there seems to be a
    // performance hit above
    maxMassConfigurationSkipping = runningMaxClusterAlgorithmMass-2;
    for(G4int i=0; i<runningMaxClusterAlgorithmMass-2; ++i) // no caching for A=1,2
      checkedConfigurations[i].clear();
#endif

    // Initialise position, momentum and energy of the running cluster
    // configuration
    runningPositions[1] = leadingParticlePosition;
    runningMomenta[1] = leadingParticleMomentum;
    runningEnergies[1] = theLeadingParticle->getEnergy();
    runningPotentials[1] = theLeadingParticle->getPotentialEnergy();

    // Make sure that all the elements of isInRunningConfiguration are false.
// assert(std::count(isInRunningConfiguration, isInRunningConfiguration+nConsidered, true)==0);

    // Start the cluster search!
    findClusterStartingFrom(1, theLeadingParticle->getZ());

    // Again, make sure that all the elements of isInRunningConfiguration have
    // been reset to false. This is a sanity check.
// assert(std::count(isInRunningConfiguration, isInRunningConfiguration+nConsidered, true)==0);

    Cluster *chosenCluster = 0;
    if(selectedA!=0) { // A cluster was found!
      candidateConfiguration[selectedA-1] = theLeadingParticle;
      chosenCluster =  new Cluster(candidateConfiguration,
          candidateConfiguration + selectedA);
    }

    // Restore the original position of the leading particle
    theLeadingParticle->setPosition(oldLeadingParticlePosition);

    return chosenCluster;
  }

  G4double ClusteringModelIntercomparison::getPhaseSpace(const G4int oldA, ConsideredPartner const &p) {
    const G4double psSpace = (p.position - runningPositions[oldA]).mag2();
    const G4double psMomentum = (p.momentum*oldA - runningMomenta[oldA]).mag2();
    return psSpace * psMomentum * clusterPosFact2[oldA + 1];
  }

  void ClusteringModelIntercomparison::findClusterStartingFrom(const G4int oldA, const G4int oldZ) {
    const G4int newA = oldA + 1;
    const G4int oldAMinusOne = oldA - 1;
    G4int newZ;
    G4int newN;

    // Look up the phase-space cut
    const G4double phaseSpaceCut = clusterPhaseSpaceCut[newA];

    // Configuration caching enabled only for a certain mass interval
    const G4bool cachingEnabled = (newA<=maxMassConfigurationSkipping && newA>=3);

    // Set the pointer to the container of cached configurations
#if defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_HashMask)
    HashContainer *theHashContainer;
    if(cachingEnabled)
      theHashContainer = &(checkedConfigurations[oldA-2]);
    else
      theHashContainer = NULL;
#elif defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_Set)
    SortedNucleonConfigurationContainer *theConfigContainer;
    if(cachingEnabled)
      theConfigContainer = &(checkedConfigurations[oldA-2]);
    else
      theConfigContainer = NULL;
#elif !defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_None)
#error Unrecognized INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON. Allowed values are: Set, HashMask, None.
#endif

    // Minimum and maximum Z values for this mass
    const G4int ZMinForNewA = clusterZMin[newA];
    const G4int ZMaxForNewA = clusterZMax[newA];

    for(G4int i=0; i<nConsidered; ++i) {
      // Only accept particles that are not already part of the cluster
      if(isInRunningConfiguration[i]) continue;

      ConsideredPartner const &candidateNucleon = consideredPartners[i];

      // Z and A of the new cluster
      newZ = oldZ + candidateNucleon.Z;
      newN = newA - newZ;

      // Skip this nucleon if we already have too many protons or neutrons
      if(newZ > clusterZMaxAll || newN > clusterNMaxAll)
        continue;

      // Compute the phase space factor for a new cluster which
      // consists of the previous running cluster and the new
      // candidate nucleon:
      const G4double phaseSpace = getPhaseSpace(oldA, candidateNucleon);
      if(phaseSpace > phaseSpaceCut) continue;

      // Store the candidate nucleon in the running configuration
      runningConfiguration[oldAMinusOne] = i;
#if defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_HashMask)
      Hashing::HashType configHash;
      HashIterator aHashIter;
#elif defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_Set)
      SortedNucleonConfiguration thisConfig;
      SortedNucleonConfigurationIterator thisConfigIter;
#endif
      if(cachingEnabled) {
#if defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_HashMask)
        configHash = Hashing::hashConfig(runningConfiguration, oldA);
        aHashIter = theHashContainer->lower_bound(configHash);
        // If we have already checked this configuration, skip it
        if(aHashIter!=theHashContainer->end()
           && !(configHash < *aHashIter))
          continue;
#elif defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_Set)
        thisConfig.fill(runningConfiguration,oldA);
        thisConfigIter = theConfigContainer->lower_bound(thisConfig);
        // If we have already checked this configuration, skip it
        if(thisConfigIter!=theConfigContainer->end()
           && !(thisConfig < *thisConfigIter))
          continue;
#endif
      }

      // Sum of the total energies of the cluster components
      runningEnergies[newA] = runningEnergies[oldA] + candidateNucleon.energy;
      // Sum of the potential energies of the cluster components
      runningPotentials[newA] = runningPotentials[oldA] + candidateNucleon.potentialEnergy;

      // Update the available cascading kinetic energy
      G4double oldCascadingEnergyPool = cascadingEnergyPool;
      if(!candidateNucleon.isTargetSpectator)
        cascadingEnergyPool -= candidateNucleon.energy - candidateNucleon.potentialEnergy - 931.3;

      // Check an approximate Coulomb barrier. If the cluster is below
      // 0.5*barrier and the remaining available energy from cascading nucleons
      // will not bring it above, reject the cluster.
      const G4double halfB = 0.72 * newZ *
        theNucleus->getZ()/(theNucleus->getDensity()->getProtonNuclearRadius()+1.7);
      const G4double tout = runningEnergies[newA] - runningPotentials[newA] -
        931.3*newA;
      if(tout<=halfB && tout+cascadingEnergyPool<=halfB) {
        cascadingEnergyPool = oldCascadingEnergyPool;
        continue;
      }

      // Here the nucleon has passed all the tests. Accept it in the cluster.
      runningPositions[newA] = (runningPositions[oldA] * oldA + candidateNucleon.position)*clusterPosFact[newA];
      runningMomenta[newA] = runningMomenta[oldA] + candidateNucleon.momentum;

      // Add the config to the container
      if(cachingEnabled) {
#if defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_HashMask)
        theHashContainer->insert(aHashIter, configHash);
#elif defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_Set)
        theConfigContainer->insert(thisConfigIter, thisConfig);
#endif
      }

      // Set the flag that reminds us that this nucleon has already been taken
      // in the running configuration
      isInRunningConfiguration[i] = true;

      // Keep track of the best physical cluster
      if(newZ >= ZMinForNewA && newZ <= ZMaxForNewA) {
        // Note: sqc is real kinetic energy, not the square of the kinetic energy!
        const G4double sqc = KinematicsUtils::invariantMass(runningEnergies[newA],
            runningMomenta[newA]);
        const G4double sqct = (sqc - 2.*newZ*protonMass - 2.*(newA-newZ)*neutronMass
             + ParticleTable::getRealMass(newA, newZ))
          *clusterPosFact[newA];

        if(sqct < sqtot) {
          // This is the best cluster we have found so far. Store its
          // kinematics.
          sqtot = sqct;
          selectedA = newA;
          selectedZ = newZ;

          // Store the running configuration in a ParticleList
          for(G4int j=0; j<oldA; ++j)
            candidateConfiguration[j] = consideredPartners[runningConfiguration[j]].particle;

          // Sanity check on number of nucleons in running configuration
// assert(std::count(isInRunningConfiguration, isInRunningConfiguration+nConsidered, true)==selectedA-1);
        }
      }

      // The method recursively calls itself for the next mass
      if(newA < runningMaxClusterAlgorithmMass && newA+1 < theNucleus->getA()) {
	findClusterStartingFrom(newA, newZ);
      }

      // Reset the running configuration flag and the cascading energy pool
      isInRunningConfiguration[i] = false;
      cascadingEnergyPool = oldCascadingEnergyPool;
    }
  }

  G4bool ClusteringModelIntercomparison::clusterCanEscape(Nucleus const * const n, Cluster const * const c) {
    // Forbid emission of the whole nucleus
    if(c->getA()>=n->getA())
      return false;

    // Check the escape angle of the cluster
    const ThreeVector &pos = c->getPosition();
    const ThreeVector &mom = c->getMomentum();
    const G4double cosEscapeAngle = pos.dot(mom) / std::sqrt(pos.mag2()*mom.mag2());
    if(cosEscapeAngle < limitCosEscapeAngle)
      return false;

    return true;
  }

}
