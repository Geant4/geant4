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

#include "G4INCLReflectionChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLRandom.hh"
#include "G4INCLINuclearPotential.hh"

#include <cmath>

namespace G4INCL {
  ReflectionChannel::ReflectionChannel(Nucleus *n, Particle *p)
    :theNucleus(n),theParticle(p)
  {
  }

  ReflectionChannel::~ReflectionChannel()
  {
  }

  FinalState* ReflectionChannel::getFinalState()
  {
    FinalState *fs = new FinalState(); // Create final state for the output
    fs->setTotalEnergyBeforeInteraction(theParticle->getEnergy() - theParticle->getPotentialEnergy());

    G4double pspr = theParticle->getPosition().dot(theParticle->getMomentum());
    G4double x2cour = theParticle->getPosition().mag2();
    ThreeVector newMomentum = theParticle->getMomentum() - (theParticle->getPosition() * (2.0 * pspr/x2cour));
    //ThreeVector newMomentum = -theParticle->getMomentum(); // For debugging
    theParticle->setMomentum(newMomentum);
    theNucleus->updatePotentialEnergy(theParticle);
    fs->addModifiedParticle(theParticle);

    return fs;
  }
}

