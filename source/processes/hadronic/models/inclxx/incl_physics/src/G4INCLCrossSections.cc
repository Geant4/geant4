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

#include "G4INCLCrossSections.hh"
#include "G4INCLKinematicsUtils.hh"
#include "G4INCLParticleTable.hh"
#include "G4INCLLogger.hh"
#include "G4INCLCrossSectionsINCL46.hh"
#include "G4INCLCrossSectionsMultiPions.hh"
// #include <cassert>

namespace G4INCL {

  namespace {
    G4ThreadLocal ICrossSections *theCrossSections;
  }

  namespace CrossSections {
    G4double elastic(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->elastic(p1,p2);
    }

    G4double total(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->total(p1,p2);
    }

    G4double NDeltaToNN(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NDeltaToNN(p1,p2);
    }

    G4double NNToNDelta(Particle const * const p1, Particle const * const p2) {
      return theCrossSections->NNToNDelta(p1,p2);
    }

      G4double NNToxPiNN(const G4int xpi, Particle const * const p1, Particle const * const p2) {
          return theCrossSections->NNToxPiNN(xpi,p1,p2);
      }

      G4double piNToDelta(Particle const * const p1, Particle const * const p2) {
          return theCrossSections->piNToDelta(p1,p2);
      }

      G4double piNToxPiN(const G4int xpi, Particle const * const p1, Particle const * const p2) {
          return theCrossSections->piNToxPiN(xpi,p1,p2);
      }

    G4double calculateNNAngularSlope(G4double energyCM, G4int iso) {
      return theCrossSections->calculateNNAngularSlope(energyCM, iso);
    }

    G4double interactionDistancePiN(const G4double projectileKineticEnergy) {
      ThreeVector nullVector;
      ThreeVector unitVector(0., 0., 1.);

      Particle piPlusProjectile(PiPlus, unitVector, nullVector);
      piPlusProjectile.setEnergy(piPlusProjectile.getMass()+projectileKineticEnergy);
      piPlusProjectile.adjustMomentumFromEnergy();
      Particle piZeroProjectile(PiZero, unitVector, nullVector);
      piZeroProjectile.setEnergy(piZeroProjectile.getMass()+projectileKineticEnergy);
      piZeroProjectile.adjustMomentumFromEnergy();
      Particle piMinusProjectile(PiMinus, unitVector, nullVector);
      piMinusProjectile.setEnergy(piMinusProjectile.getMass()+projectileKineticEnergy);
      piMinusProjectile.adjustMomentumFromEnergy();

      Particle protonTarget(Proton, nullVector, nullVector);
      Particle neutronTarget(Neutron, nullVector, nullVector);
      const G4double sigmapipp = total(&piPlusProjectile, &protonTarget);
      const G4double sigmapipn = total(&piPlusProjectile, &neutronTarget);
      const G4double sigmapi0p = total(&piZeroProjectile, &protonTarget);
      const G4double sigmapi0n = total(&piZeroProjectile, &neutronTarget);
      const G4double sigmapimp = total(&piMinusProjectile, &protonTarget);
      const G4double sigmapimn = total(&piMinusProjectile, &neutronTarget);
      /* We compute the interaction distance from the largest of the pi-N cross
       * sections. Note that this is different from INCL4.6, which just takes the
       * average of the six, and will in general lead to a different geometrical
       * cross section.
       */
      const G4double largestSigma = std::max(sigmapipp, std::max(sigmapipn, std::max(sigmapi0p, std::max(sigmapi0n, std::max(sigmapimp,sigmapimn)))));
      const G4double interactionDistance = std::sqrt(largestSigma/Math::tenPi);

      return interactionDistance;
    }

    G4double interactionDistanceNN(const ParticleSpecies &aSpecies, const G4double kineticEnergy) {
// assert(aSpecies.theType==Proton || aSpecies.theType==Neutron || aSpecies.theType==Composite);
// assert(aSpecies.theA>0);
      ThreeVector nullVector;
      ThreeVector unitVector(0.,0.,1.);

      const G4double kineticEnergyPerNucleon = kineticEnergy / aSpecies.theA;

      Particle protonProjectile(Proton, unitVector, nullVector);
      protonProjectile.setEnergy(protonProjectile.getMass()+kineticEnergyPerNucleon);
      protonProjectile.adjustMomentumFromEnergy();
      Particle neutronProjectile(Neutron, unitVector, nullVector);
      neutronProjectile.setEnergy(neutronProjectile.getMass()+kineticEnergyPerNucleon);
      neutronProjectile.adjustMomentumFromEnergy();

      Particle protonTarget(Proton, nullVector, nullVector);
      Particle neutronTarget(Neutron, nullVector, nullVector);
      const G4double sigmapp = total(&protonProjectile, &protonTarget);
      const G4double sigmapn = total(&protonProjectile, &neutronTarget);
      const G4double sigmann = total(&neutronProjectile, &neutronTarget);
      /* We compute the interaction distance from the largest of the NN cross
       * sections. Note that this is different from INCL4.6, which just takes the
       * average of the four, and will in general lead to a different geometrical
       * cross section.
       */
      const G4double largestSigma = std::max(sigmapp, std::max(sigmapn, sigmann));
      const G4double interactionDistance = std::sqrt(largestSigma/Math::tenPi);

      return interactionDistance;
    }

    void setCrossSections(ICrossSections *c) {
      theCrossSections = c;
    }

    void deleteCrossSections() {
      delete theCrossSections;
      theCrossSections = NULL;
    }

    void initialize(Config const * const theConfig) {
      CrossSectionsType crossSections = theConfig->getCrossSectionsType();
      if(crossSections == INCL46CrossSections)
        setCrossSections(new CrossSectionsINCL46);
      else if(crossSections == MultiPionsCrossSections)
        setCrossSections(new CrossSectionsMultiPions);
    }
  }
}
