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

/*
 * Particle.cc
 *
 *  \date Jun 5, 2009
 * \author Pekka Kaitaniemi
 */

#include "G4INCLParticle.hh"
#include "G4INCLParticleTable.hh"

namespace G4INCL {

  long Particle::nextID = 1;

  Particle::Particle()
    : theZ(0), theA(0),
    theParticipantType(TargetSpectator),
    theType(UnknownParticle),
    theEnergy(0.0),
    thePropagationEnergy(&theEnergy),
    theFrozenEnergy(theEnergy),
    theMomentum(ThreeVector(0.,0.,0.)),
    thePropagationMomentum(&theMomentum),
    theFrozenMomentum(theMomentum),
    thePosition(ThreeVector(0.,0.,0.)),
    nCollisions(0),
    nDecays(0),
    thePotentialEnergy(0.0),
    theHelicity(0.0),
    emissionTime(0.0),
    outOfWell(false),
    theMass(0.)
  {
    ID = nextID;
    nextID++;
  }

  Particle::Particle(ParticleType t, G4double energy,
      ThreeVector momentum, ThreeVector position)
    : theEnergy(energy),
    thePropagationEnergy(&theEnergy),
    theFrozenEnergy(theEnergy),
    theMomentum(momentum),
    thePropagationMomentum(&theMomentum),
    theFrozenMomentum(theMomentum),
    thePosition(position),
    nCollisions(0), nDecays(0),
      thePotentialEnergy(0.), theHelicity(0.0),
      emissionTime(0.0), outOfWell(false)
  {
    theParticipantType = TargetSpectator;
    ID = nextID;
    nextID++;
    if(theEnergy <= 0.0) {
      WARN("Particle with energy " << theEnergy << " created." << std::endl);
    }
    setType(t);
    setMass(getInvariantMass());
  }

  Particle::Particle(ParticleType t,
      ThreeVector momentum, ThreeVector position)
    : thePropagationEnergy(&theEnergy),
    theMomentum(momentum),
    thePropagationMomentum(&theMomentum),
    theFrozenMomentum(theMomentum),
    thePosition(position),
    nCollisions(0), nDecays(0),
      thePotentialEnergy(0.), theHelicity(0.0),
      emissionTime(0.0), outOfWell(false)
  {
    theParticipantType = TargetSpectator;
    ID = nextID;
    nextID++;
    setType(t);
    if( isResonance() ) {
      ERROR("Cannot create resonance without specifying its momentum four-vector." << std::endl);
    }
    G4double energy = std::sqrt(theMomentum.mag2() + theMass*theMass);
    theEnergy = energy;
    theFrozenEnergy = theEnergy;
  }

  const ThreeVector &Particle::adjustMomentumFromEnergy() {
    const G4double p2 = theMomentum.mag2();
    G4double newp2 = theEnergy*theEnergy - theMass*theMass;
    if( newp2<0.0 ) {
      ERROR("Particle has E^2 < m^2." << std::endl << print());
      newp2 = 0.0;
      theEnergy = theMass;
    }

    theMomentum *= std::sqrt(newp2/p2);
    return theMomentum;
  }

  G4double Particle::adjustEnergyFromMomentum() {
    theEnergy = std::sqrt(theMomentum.mag2() + theMass*theMass);
    return theEnergy;
  }
}
