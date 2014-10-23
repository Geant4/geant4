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

#ifndef G4INCLClusteringModelIntercomparison_hh
#define G4INCLClusteringModelIntercomparison_hh 1

#ifdef INCLXX_IN_GEANT4_MODE
#define INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_Set 1
#endif // INCLXX_IN_GEANT4_MODE

#include "G4INCLIClusteringModel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLCluster.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLHashing.hh"

#include <set>
#include <algorithm>

namespace G4INCL {

  /** \brief Container for the relevant information
   *
   * This struct contains all the information that is relevant for the
   * clustering algorithm. It is probably more compact than the Particles it
   * feeds on, hopefully improving cache performance.
   */
  struct ConsideredPartner {
    Particle *particle;
    G4bool isTargetSpectator;
    G4int Z;
    ThreeVector position;
    ThreeVector momentum;
    G4double energy;
    G4double potentialEnergy;

    ConsideredPartner() :
      particle(NULL),
      isTargetSpectator(false),
      Z(0),
      energy(0.),
      potentialEnergy(0.)
    {}

    ConsideredPartner(Particle * const p) :
      particle(p),
      isTargetSpectator(particle->isTargetSpectator()),
      Z(particle->getZ()),
      position(particle->getPosition()),
      momentum(particle->getMomentum()),
      energy(particle->getEnergy()),
      potentialEnergy(particle->getPotentialEnergy())
    {}
  };

  /// \brief Cluster coalescence algorithm used in the IAEA intercomparison
  class ClusteringModelIntercomparison : public IClusteringModel {
  public:
    ClusteringModelIntercomparison(Config const * const theConfig) :
      theNucleus(NULL),
      selectedA(0),
      selectedZ(0),
      sqtot(0.),
      cascadingEnergyPool(0.),
      protonMass(ParticleTable::getRealMass(Proton)),
      neutronMass(ParticleTable::getRealMass(Neutron)),
      runningMaxClusterAlgorithmMass(theConfig->getClusterMaxMass()),
      nConsideredMax(0),
      nConsidered(0),
      consideredPartners(NULL),
      isInRunningConfiguration(NULL),
      maxMassConfigurationSkipping(ParticleTable::maxClusterMass)
    {
      // Set up the maximum charge and neutron number for clusters
      clusterZMaxAll = 0;
      clusterNMaxAll = 0;
      for(G4int A=0; A<=runningMaxClusterAlgorithmMass; ++A) {
        if(clusterZMax[A]>clusterZMaxAll)
          clusterZMaxAll = clusterZMax[A];
        if(A-clusterZMin[A]>clusterNMaxAll)
          clusterNMaxAll = A-clusterZMin[A];
      }
      std::fill(candidateConfiguration,
                candidateConfiguration + ParticleTable::maxClusterMass,
                static_cast<Particle*>(NULL));

      std::fill(runningEnergies,
                runningEnergies + ParticleTable::maxClusterMass,
                0.0);

      std::fill(runningPotentials,
                runningPotentials + ParticleTable::maxClusterMass,
                0.0);

      std::fill(runningConfiguration,
                runningConfiguration + ParticleTable::maxClusterMass,
                -1);

    }

    virtual ~ClusteringModelIntercomparison() {
      delete [] consideredPartners;
      delete [] isInRunningConfiguration;
    }

    virtual Cluster* getCluster(Nucleus*, Particle*);
    virtual G4bool clusterCanEscape(Nucleus const * const, Cluster const * const);

  private:
    void findClusterStartingFrom(const G4int oldA, const G4int oldZ);
    G4double getPhaseSpace(const G4int oldA, ConsideredPartner const &p);

    Nucleus *theNucleus;

    G4double runningEnergies[ParticleTable::maxClusterMass+1];
    ThreeVector runningMomenta[ParticleTable::maxClusterMass+1];
    ThreeVector runningPositions[ParticleTable::maxClusterMass+1];
    G4double runningPotentials[ParticleTable::maxClusterMass+1];
#if defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_HashMask)
    Hashing::NucleonItem runningConfiguration[ParticleTable::maxClusterMass];
#elif defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_Set) || defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_None)
    G4int runningConfiguration[ParticleTable::maxClusterMass];
#else
#error Unrecognized INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON. Allowed values are: Set, HashMask, None.
#endif

    G4int selectedA, selectedZ;
    G4double sqtot;

    G4int clusterZMaxAll, clusterNMaxAll;

    G4double cascadingEnergyPool;

    /// \brief Lower limit of Z for cluster of mass A
    static const G4int clusterZMin[ParticleTable::maxClusterMass+1];
    /// \brief Upper limit of Z for cluster of mass A
    static const G4int clusterZMax[ParticleTable::maxClusterMass+1];

    /// \brief Precomputed factor 1.0/A
    static const G4double clusterPosFact[ParticleTable::maxClusterMass+1];

    /// \brief Precomputed factor (1.0/A)^2
    static const G4double clusterPosFact2[ParticleTable::maxClusterMass+1];

    /// \brief Phase-space parameters for cluster formation
    static const G4double clusterPhaseSpaceCut[ParticleTable::maxClusterMass+1];

    static const G4double limitCosEscapeAngle;

    const G4double protonMass;
    const G4double neutronMass;

    G4int runningMaxClusterAlgorithmMass;

    G4int nConsideredMax;
    G4int nConsidered;

    /** \brief Array of considered cluster partners
     *
     * A dynamical array of ConsideredPartner objects is allocated on this
     * variable and filled with pointers to nucleons which are eligible for
     * clustering. We used to use a ParticleList for this purpose, but this
     * made it very cumbersome to check whether nucleons had already been
     * included in the running configuration. Using an array of Particle*
     * coupled with a boolean mask (\see{isInRunningConfiguration}) reduces the
     * overhead by a large amount.  Running times for 1-GeV p+Pb208 went down
     * by almost 30% (!).
     *
     * Lesson learnt: when you need speed, nothing beats a good ol' array.
     */
    ConsideredPartner *consideredPartners;

    /** \brief Array of flags for nucleons in the running configuration
     *
     * Clustering partners that are already used in the running cluster
     * configuration are flagged as "true" in this array.
     */
    G4bool *isInRunningConfiguration;

    /** \brief Best cluster configuration
     *
     * This array contains pointers to the nucleons which make up the best
     * cluster configuration that has been found so far.
     */
    Particle *candidateConfiguration[ParticleTable::maxClusterMass];

#if defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_HashMask)
    typedef std::set<Hashing::HashType> HashContainer;
    typedef HashContainer::iterator HashIterator;

    /// \brief Array of containers for configurations that have already been checked
    HashContainer checkedConfigurations[ParticleTable::maxClusterMass-2];
#elif defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_Set)
    /** \brief Class for storing and comparing sorted nucleon configurations
     *
     * This class is actually just a wrapper around an array of Particle*
     * pointers. It provides a lexicographical comparison operator
     * (SortedNucleonConfiguration::operator<) for inclusion in std::set
     * containers.
     */
    class SortedNucleonConfiguration {
      public:
        // Use Particle* as nucleon identifiers
        typedef G4int NucleonItem;

        /// \brief Constructor
        SortedNucleonConfiguration() : theSize(0), nucleons(NULL) {}

        /// \brief Copy constructor
        SortedNucleonConfiguration(const SortedNucleonConfiguration &rhs) :
          theSize(rhs.theSize),
          nucleons(new NucleonItem[theSize])
      {
        std::copy(rhs.nucleons, rhs.nucleons+theSize, nucleons);
      }

        /// \brief Destructor
        ~SortedNucleonConfiguration() {
          delete [] nucleons;
        }

        /// \brief Helper method for the assignment operator
        void swap(SortedNucleonConfiguration &rhs) {
          std::swap(theSize, rhs.theSize);
          std::swap(nucleons, rhs.nucleons);
        }

        /// \brief Assignment operator
        SortedNucleonConfiguration &operator=(const SortedNucleonConfiguration &rhs) {
          SortedNucleonConfiguration tempConfig(rhs);
          swap(tempConfig);
          return *this;
        }

        /** \brief Order operator for SortedNucleonConfiguration
         *
         * The comparison is done lexicographically (i.e. from the first
         * element to the last).
         */
        G4bool operator<(const SortedNucleonConfiguration &rhs) const {
// assert(theSize==rhs.theSize);
          return std::lexicographical_compare(nucleons, nucleons+theSize, rhs.nucleons, rhs.nucleons+theSize);
        }

        /// \brief Fill configuration with array of NucleonItem
        void fill(NucleonItem *config, size_t n) {
          theSize = n;
          nucleons = new NucleonItem[theSize];
          std::copy(config, config+theSize, nucleons);
          std::sort(nucleons, nucleons+theSize);
        }

      private:
        /// \brief Size of the array
        size_t theSize;

        /// \brief The real array
        NucleonItem *nucleons;
    };

    typedef std::set<SortedNucleonConfiguration> SortedNucleonConfigurationContainer;
    typedef SortedNucleonConfigurationContainer::iterator SortedNucleonConfigurationIterator;

    /// \brief Array of containers for configurations that have already been checked
    SortedNucleonConfigurationContainer checkedConfigurations[ParticleTable::maxClusterMass-2];
#elif !defined(INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON_None)
#error Unrecognized INCL_CACHING_CLUSTERING_MODEL_INTERCOMPARISON. Allowed values are: Set, HashMask, None.
#endif

    /** \brief Maximum mass for configuration storage
     *
     * Skipping configurations becomes inefficient above this mass.
     */
    G4int maxMassConfigurationSkipping;
  };

}

#endif
