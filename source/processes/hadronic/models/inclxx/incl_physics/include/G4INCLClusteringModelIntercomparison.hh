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
// INCL++ revision: v5.1.5
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#ifndef G4INCLClusteringModelIntercomparison_hh
#define G4INCLClusteringModelIntercomparison_hh 1

#include "G4INCLIClusteringModel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLCluster.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLKinematicsUtils.hh"

namespace G4INCL {

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
      isInRunningConfiguration(NULL)
    {
      // Set up the maximum charge and neutron number for clusters
      clusterZMaxAll = 0;
      clusterNMaxAll = 0;
      for(G4int A=0; A<=runningMaxClusterAlgorithmMass; ++A) {
        if(ParticleTable::clusterZMax[A]>clusterZMaxAll)
          clusterZMaxAll = ParticleTable::clusterZMax[A];
        if(A-ParticleTable::clusterZMin[A]>clusterNMaxAll)
          clusterNMaxAll = A-ParticleTable::clusterZMin[A];
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
    }

    virtual ~ClusteringModelIntercomparison() {
      delete [] consideredPartners;
      delete [] isInRunningConfiguration;
    }

    virtual Cluster* getCluster(Nucleus*, Particle*);
    virtual G4bool clusterCanEscape(Nucleus const * const, Cluster const * const);

  private:
    void findClusterStartingFrom(const G4int oldA, const G4int oldZ);
    G4double getPhaseSpace(const G4int oldA, Particle const * const p);

    Nucleus *theNucleus;

    G4double runningEnergies[ParticleTable::maxClusterMass];
    ThreeVector runningMomenta[ParticleTable::maxClusterMass];
    ThreeVector runningPositions[ParticleTable::maxClusterMass];
    G4double runningPotentials[ParticleTable::maxClusterMass];

    G4int selectedA, selectedZ;
    G4double sqtot;

    G4int clusterZMaxAll, clusterNMaxAll;

    G4double cascadingEnergyPool;

    static const G4double limitCosEscapeAngle;

    const G4double protonMass;
    const G4double neutronMass;

    G4int runningMaxClusterAlgorithmMass;

    G4int nConsideredMax;
    G4int nConsidered;

    /** \brief Array of considered cluster partners
     *
     * A dynamical array of Particle* is allocated on this variable and filled
     * with pointers to nucleons which are eligible for clustering. We used to
     * use a ParticleList for this purpose, but this made it very cumbersome to
     * check whether nucleons had already been included in the running
     * configuration. Using an array of Particle* coupled with a boolean mask
     * (\see{isInRunningConfiguration}) reduces the overhead by a large amount.
     * Running times for 1-GeV p+Pb208 went down by almost 30% (!).
     *
     * Lesson learnt: when you need speed, nothing beats a good ol' array.
     */
    Particle **consideredPartners;

    /** \brief Array of flags for nucleons in the running configuration
     *
     * Clustering partners that are already used in the running cluster
     * configuration are flagged as "true" in this array.
     */
    G4bool *isInRunningConfiguration;

    /** \brief Best cluster configuration
     *
     * A dynamical array of Particle* is allocated on this variable and filled
     * with pointers to the nucleons which make up the best cluster
     * configuration that has been found so far.
     */
    Particle *candidateConfiguration[ParticleTable::maxClusterMass];
  };

}

#endif
