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

#ifndef G4INCLClusteringModelIntercomparison_hh
#define G4INCLClusteringModelIntercomparison_hh 1

#include "G4INCLIClusteringModel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLCluster.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLKinematicsUtils.hh"

namespace G4INCL {

  class ClusteringModelIntercomparison : public IClusteringModel {
  public:
    ClusteringModelIntercomparison() {
      zeroOut();

      // Set up the maximum charge and neutron number for clusters
      clusterZMaxAll = 0;
      clusterNMaxAll = 0;
      for(G4int A=0; A<=maxClusterAlgorithmMass; ++A) {
        if(ParticleTable::clusterZMax[A]>clusterZMaxAll)
          clusterZMaxAll = ParticleTable::clusterZMax[A];
        if(A-ParticleTable::clusterZMin[A]>clusterNMaxAll)
          clusterNMaxAll = A-ParticleTable::clusterZMin[A];
      }
    }

    void cleanUp() {
      delete candidateConfiguration;
      consideredPartners.clear();
      runningConfiguration.clear();
    }

    void zeroOut() {
      candidateConfiguration = 0;
    }

    virtual ~ClusteringModelIntercomparison() {
      cleanUp();
    }

    virtual Cluster* getCluster(Nucleus*, Particle*);
    virtual G4bool clusterCanEscape(Cluster const * const);

  private:
    void findClusterStartingFrom(const G4int oldA, const G4int oldZ);
    G4double getPhaseSpace(G4int oldA, Particle *p);

    Nucleus *theNucleus;
    Particle *theLeadingParticle;
    ParticleList consideredPartners;
    ParticleList* candidateConfiguration;

    G4double runningEnergies[ParticleTable::maxClusterMass];
    ThreeVector runningMomenta[ParticleTable::maxClusterMass];
    ThreeVector runningPositions[ParticleTable::maxClusterMass];
    ParticleList runningConfiguration; // Use deque instead?
    G4double runningPotentials[ParticleTable::maxClusterMass];

    G4int selectedA, selectedZ;
    G4double sqtot;

    G4int clusterZMaxAll, clusterNMaxAll;

    G4double participantEnergyPool;

    static const G4double limitCosEscapeAngle;
  };

}

#endif
