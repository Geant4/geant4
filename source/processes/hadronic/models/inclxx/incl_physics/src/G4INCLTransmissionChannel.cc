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

#include "G4INCLTransmissionChannel.hh"

namespace G4INCL {

  TransmissionChannel::TransmissionChannel(Nucleus *nucleus, Particle *particle)
    : theNucleus(nucleus), theParticle(particle)
  {}

  TransmissionChannel::~TransmissionChannel() {}

  G4double TransmissionChannel::particleLeaves() {

    // \todo{this is the place to add refraction}

    // The particle energy outside the nucleus. Subtract the nuclear
    // potential from the kinetic energy when leaving the nucleus
    G4double kineticEnergyOutside = theParticle->getEnergy()
      - theParticle->getPotentialEnergy()
      - theParticle->getMass();

    // Correction for real masses
    const G4int AParent = theNucleus->getA();
    const G4int ZParent = theNucleus->getZ();
    const G4double theQValueCorrection = theParticle->getEmissionQValueCorrection(AParent,ZParent);
    kineticEnergyOutside += theQValueCorrection;
    theParticle->setTableMass();

    // Scaling factor for the particle momentum
    theParticle->setEnergy(kineticEnergyOutside + theParticle->getMass());
    theParticle->adjustMomentumFromEnergy();
    theParticle->setPotentialEnergy(0.);

    return theQValueCorrection;
  }

  FinalState* TransmissionChannel::getFinalState() {
    FinalState *fs = new FinalState;
    G4double initialEnergy = 0.0;
    initialEnergy = theParticle->getEnergy() - theParticle->getPotentialEnergy();

    // Correction for real masses
    initialEnergy += theParticle->getTableMass() - theParticle->getMass();

    const G4double theQValueCorrection = particleLeaves();
    fs->setTotalEnergyBeforeInteraction(initialEnergy + theQValueCorrection);
    fs->addOutgoingParticle(theParticle); // We write the particle down as outgoing
    return fs;
  }
}
