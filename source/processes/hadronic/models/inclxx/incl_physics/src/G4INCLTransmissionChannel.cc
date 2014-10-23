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

#include "G4INCLTransmissionChannel.hh"

namespace G4INCL {

  TransmissionChannel::TransmissionChannel(Nucleus * const nucleus, Particle * const particle)
    : theNucleus(nucleus), theParticle(particle),
    refraction(false),
    pOutMag(0.),
    kineticEnergyOutside(initializeKineticEnergyOutside()),
    cosRefractionAngle(1.)
  {}

  TransmissionChannel::TransmissionChannel(Nucleus * const nucleus, Particle * const particle, const G4double TOut)
    : theNucleus(nucleus), theParticle(particle),
    refraction(false),
    pOutMag(0.),
    kineticEnergyOutside(TOut),
    cosRefractionAngle(1.)
  {}

  TransmissionChannel::TransmissionChannel(Nucleus * const nucleus, Particle * const particle, const G4double kOut, const G4double cosR)
    : theNucleus(nucleus), theParticle(particle),
    refraction(true),
    pOutMag(kOut),
    kineticEnergyOutside(initializeKineticEnergyOutside()),
    cosRefractionAngle(cosR)
  {}

  TransmissionChannel::~TransmissionChannel() {}

  G4double TransmissionChannel::initializeKineticEnergyOutside() {
    // The particle energy outside the nucleus. Subtract the nuclear
    // potential from the kinetic energy when leaving the nucleus
    G4double TOut = theParticle->getEnergy()
      - theParticle->getPotentialEnergy()
      - theParticle->getMass();

    // Correction for real masses
    const G4int AParent = theNucleus->getA();
    const G4int ZParent = theNucleus->getZ();
    const G4double theQValueCorrection = theParticle->getEmissionQValueCorrection(AParent,ZParent);
    TOut += theQValueCorrection;
    return TOut;
  }

  void TransmissionChannel::particleLeaves() {

    // Use the table mass in the outside world
    theParticle->setTableMass();
    theParticle->setPotentialEnergy(0.);

    if(refraction) {
      // Change the momentum direction
      // The magnitude of the particle momentum outside the nucleus will be
      // fixed by the kineticEnergyOutside variable. This is done in order to
      // avoid numerical inaccuracies.
      const ThreeVector &position = theParticle->getPosition();
      const G4double r2 = position.mag2();
      ThreeVector normal;
      if(r2>0.)
        normal = position / std::sqrt(r2);

      const ThreeVector &momentum = theParticle->getMomentum();

      const ThreeVector pOut = normal * (pOutMag * cosRefractionAngle) + momentum - normal * normal.dot(momentum);
// assert(std::fabs(pOut.mag()-pOutMag)<1.e-5);

      theParticle->setMomentum(pOut);
    }
    // Scaling factor for the particle momentum
    theParticle->setEnergy(kineticEnergyOutside + theParticle->getMass());
    theParticle->adjustMomentumFromEnergy();
  }

  void TransmissionChannel::fillFinalState(FinalState *fs) {
    G4double initialEnergy = 0.0;
    initialEnergy = theParticle->getEnergy() - theParticle->getPotentialEnergy();

    // Correction for real masses
    const G4int AParent = theNucleus->getA();
    const G4int ZParent = theNucleus->getZ();
    initialEnergy += theParticle->getTableMass() - theParticle->getMass()
      + theParticle->getEmissionQValueCorrection(AParent,ZParent);

    particleLeaves();

    fs->setTotalEnergyBeforeInteraction(initialEnergy);
    fs->addOutgoingParticle(theParticle); // We write the particle down as outgoing
  }
}
