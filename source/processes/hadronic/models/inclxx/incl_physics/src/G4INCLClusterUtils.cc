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

#include "G4INCLClusterUtils.hh"

namespace G4INCL {

  G4double ClusterUtils::getTotalEnergy(const ParticleList &pl) {
    G4double E = 0.0;
    for(ParticleIter i = pl.begin(); i != pl.end(); ++i) {
      E += (*i)->getEnergy();
    }
    return E;
  }

  G4double ClusterUtils::getKineticEnergy(const ParticleList &pl) {
    G4double Ekin = 0.0;
    for(ParticleIter i = pl.begin(); i != pl.end(); ++i) {
      Ekin += std::sqrt(std::pow((*i)->getEnergy(), 2) - std::pow((*i)->getMass(), 2));;
    }
    return Ekin;
  }

  G4int ClusterUtils::getZ(const ParticleList &pl) {
    G4int Z = 0;
    for(ParticleIter i = pl.begin(); i != pl.end(); ++i) {
      Z += (*i)->getZ();
    }
    return Z;
  }

  G4int ClusterUtils::getZ(const ParticleList &pl, Particle *p) {
    return (ClusterUtils::getZ(pl) + p->getZ());
  }

  G4int ClusterUtils::getA(const ParticleList &pl) {
    G4int A = 0;
    for(ParticleIter i = pl.begin(); i != pl.end(); ++i) {
      A += (*i)->getA();
    }
    return A;
  }

  G4int ClusterUtils::getA(const ParticleList &pl, Particle *p) {
    return (ClusterUtils::getA(pl) + p->getA());
  }

  ThreeVector ClusterUtils::getNewPositionVector(const ParticleList &pl)
  {
    ThreeVector pos(0.0, 0.0, 0.0);
    G4int A = 1;
    for(ParticleIter i = pl.begin(); i != pl.end(); ++i) {
      pos += (pos * A)*ParticleTable::clusterPosFact[A];
      ++A;
    }
    return pos;
  }

  ThreeVector ClusterUtils::getNewPositionVector(const ParticleList &pl, Particle *p)
  {
    ThreeVector pos = ClusterUtils::getNewPositionVector(pl);
    return ((pos * pl.size()) + p->getPosition()) * ParticleTable::clusterPosFact[pl.size() + 1];
  }

  ThreeVector ClusterUtils::getNewPositionVector(const ThreeVector &oldPosition, const ParticleList &pl, Particle *p) {
    ThreeVector newPosition = oldPosition;
    newPosition *= pl.size();
    newPosition += p->getPosition();
    newPosition *= ParticleTable::clusterPosFact[pl.size() + 1];
    return newPosition;
  }

  G4double ClusterUtils::getPhaseSpace(const ThreeVector &clusterPosition,
				     const ThreeVector &clusterMomentum,
				     G4int clusterA,
				     Particle *p) {
    G4double psSpace = (p->getPosition() - clusterPosition).mag2();
    G4double psMomentum = (p->getMomentum() - clusterMomentum).mag2();
    return psSpace * psMomentum * ParticleTable::clusterPosFact2[clusterA + p->getA()];
  }

  G4bool ClusterUtils::isBetterCluster(ParticleList *, ParticleList *) {
    return true;
  }
}
