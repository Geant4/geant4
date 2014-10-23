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

#include "G4INCLParticleEntryChannel.hh"
#include "G4INCLRootFinder.hh"
#include "G4INCLIntersection.hh"
#include <algorithm>

namespace G4INCL {

  ParticleEntryChannel::ParticleEntryChannel(Nucleus *n, Particle *p)
    :theNucleus(n), theParticle(p)
  {}

  ParticleEntryChannel::~ParticleEntryChannel()
  {}

  void ParticleEntryChannel::fillFinalState(FinalState *fs) {
    // Behaves slightly differency if a third body (the projectile) is present
    G4bool isNN = theNucleus->isNucleusNucleusCollision();

    /* Corrections to the energy of the entering nucleon
     *
     * In particle-nucleus reactions, the goal of this correction is to satisfy
     * energy conservation in particle-nucleus reactions using real particle
     * and nuclear masses.
     *
     * In nucleus-nucleus reactions, in addition to the above, the correction
     * is determined by a model for the excitation energy of the
     * quasi-projectile (QP). The energy of the entering nucleon is such that
     * the QP excitation energy, as determined by conservation, is what given
     * by our model.
     *
     * Possible choices for the correction (or, equivalently, for the QP
     * excitation energy):
     *
     * 1. the correction is 0. (same as in particle-nucleus);
     * 2. the correction is the separation energy of the entering nucleon in
     *    the current QP;
     * 3. the QP excitation energy is given by A. Boudard's algorithm, as
     *    implemented in INCL4.2-HI/Geant4.
     * 4. the QP excitation energy vanishes.
     *
     * Ideally, the QP excitation energy should always be >=0. Algorithms 1.
     * and 2. do not guarantee this, although violations to the rule seem to be
     * more severe for 1. than for 2.. Algorithms 3. and 4., by construction,
     * yields non-negative QP excitation energies.
     */
    G4double theCorrection;
    if(isNN) {
// assert(theParticle->isNucleon());
      ProjectileRemnant * const projectileRemnant = theNucleus->getProjectileRemnant();
// assert(projectileRemnant);

      // No correction (model 1. above)
      /*
      theCorrection = theParticle->getEmissionQValueCorrection(
          theNucleus->getA() + theParticle->getA(),
          theNucleus->getZ() + theParticle->getZ())
        + theParticle->getTableMass() - theParticle->getINCLMass();
      const G4double theProjectileCorrection = 0.;
      */

      // Correct the energy of the entering particle for the Q-value of the
      // emission from the projectile (model 2. above)
      /*
      theCorrection = theParticle->getTransferQValueCorrection(
          projectileRemnant->getA(), projectileRemnant->getZ(),
          theNucleus->getA(), theNucleus->getZ());
      G4double theProjectileCorrection;
      if(projectileRemnant->getA()>theParticle->getA()) { // if there are any particles left
        // Compute the projectile Q-value (to be used as a correction to the
        // other components of the projectile remnant)
        theProjectileCorrection = ParticleTable::getTableQValue(
            projectileRemnant->getA() - theParticle->getA(),
            projectileRemnant->getZ() - theParticle->getZ(),
            theParticle->getA(),
            theParticle->getZ());
      } else
        theProjectileCorrection = 0.;
      */

      // Fix the correction in such a way that the quasi-projectile excitation
      // energy is given by A. Boudard's INCL4.2-HI model (model 3. above).
      const G4double theProjectileExcitationEnergy =
        (projectileRemnant->getA()-theParticle->getA()>1) ?
        (projectileRemnant->computeExcitationEnergyExcept(theParticle->getID())) :
        0.;
      // Set the projectile excitation energy to zero (cold quasi-projectile,
      // model 4. above).
      // const G4double theProjectileExcitationEnergy = 0.;
      // The part that follows is common to model 3. and 4.
      const G4double theProjectileEffectiveMass =
        ParticleTable::getTableMass(projectileRemnant->getA() - theParticle->getA(), projectileRemnant->getZ() - theParticle->getZ())
        + theProjectileExcitationEnergy;
      const ThreeVector &theProjectileMomentum = projectileRemnant->getMomentum() - theParticle->getMomentum();
      const G4double theProjectileEnergy = std::sqrt(theProjectileMomentum.mag2() + theProjectileEffectiveMass*theProjectileEffectiveMass);
      const G4double theProjectileCorrection = theProjectileEnergy - (projectileRemnant->getEnergy() - theParticle->getEnergy());
      theCorrection = theParticle->getEmissionQValueCorrection(
          theNucleus->getA() + theParticle->getA(),
          theNucleus->getZ() + theParticle->getZ())
        + theParticle->getTableMass() - theParticle->getINCLMass()
        + theProjectileCorrection;
      // end of part common to model 3. and 4.


      projectileRemnant->removeParticle(theParticle, theProjectileCorrection);
    } else {
      const G4int ACN = theNucleus->getA() + theParticle->getA();
      const G4int ZCN = theNucleus->getZ() + theParticle->getZ();
      // Correction to the Q-value of the entering particle
      theCorrection = theParticle->getEmissionQValueCorrection(ACN,ZCN);
      INCL_DEBUG("The following Particle enters with correction " << theCorrection << '\n'
          << theParticle->print() << '\n');
    }

    const G4double energyBefore = theParticle->getEnergy() - theCorrection;
    G4bool success = particleEnters(theCorrection);
    fs->addEnteringParticle(theParticle);

    if(!success) {
      fs->makeParticleBelowZero();
    } else if(theParticle->isNucleon() &&
        theParticle->getKineticEnergy()<theNucleus->getPotential()->getFermiEnergy(theParticle)) {
      // If the participant is a nucleon entering below its Fermi energy, force a
      // compound nucleus
      fs->makeParticleBelowFermi();
    }

    fs->setTotalEnergyBeforeInteraction(energyBefore);
  }

  G4bool ParticleEntryChannel::particleEnters(const G4double theQValueCorrection) {

    // \todo{this is the place to add refraction}

    theParticle->setINCLMass(); // Will automatically put the particle on shell

    // Add the nuclear potential to the kinetic energy when entering the
    // nucleus

    class IncomingEFunctor : public RootFunctor {
      public:
        IncomingEFunctor(Particle * const p, Nucleus const * const n, const G4double correction) :
          RootFunctor(0., 1E6),
          theParticle(p),
          thePotential(n->getPotential()),
          theEnergy(theParticle->getEnergy()),
          theMass(theParticle->getMass()),
          theQValueCorrection(correction),
          refraction(n->getStore()->getConfig()->getRefraction()),
          theMomentumDirection(theParticle->getMomentum())
          {
            if(refraction) {
              const ThreeVector &position = theParticle->getPosition();
              const G4double r2 = position.mag2();
              if(r2>0.)
                normal = - position / std::sqrt(r2);
              G4double cosIncidenceAngle = theParticle->getCosRPAngle();
              if(cosIncidenceAngle < -1.)
                sinIncidenceAnglePOut = 0.;
              else
                sinIncidenceAnglePOut = theMomentumDirection.mag()*std::sqrt(1.-cosIncidenceAngle*cosIncidenceAngle);
            } else {
              sinIncidenceAnglePOut = 0.;
            }
          }
        ~IncomingEFunctor() {}
        G4double operator()(const G4double v) const {
          G4double energyInside = std::max(theMass, theEnergy + v - theQValueCorrection);
          theParticle->setEnergy(energyInside);
          theParticle->setPotentialEnergy(v);
          if(refraction) {
            // Compute the new direction of the particle momentum
            const G4double pIn = std::sqrt(energyInside*energyInside-theMass*theMass);
            const G4double sinRefractionAngle = sinIncidenceAnglePOut/pIn;
            const G4double cosRefractionAngle = (sinRefractionAngle>1.) ? 0. : std::sqrt(1.-sinRefractionAngle*sinRefractionAngle);
            const ThreeVector momentumInside = theMomentumDirection - normal * normal.dot(theMomentumDirection) + normal * (pIn * cosRefractionAngle);
            theParticle->setMomentum(momentumInside);
          } else {
            theParticle->setMomentum(theMomentumDirection); // keep the same direction
          }
          // Scale the particle momentum
          theParticle->adjustMomentumFromEnergy();
          return v - thePotential->computePotentialEnergy(theParticle);
        }
        void cleanUp(const G4bool /*success*/) const {}
      private:
        Particle *theParticle;
        NuclearPotential::INuclearPotential const *thePotential;
        const G4double theEnergy;
        const G4double theMass;
        const G4double theQValueCorrection;
        const G4bool refraction;
        const ThreeVector theMomentumDirection;
        ThreeVector normal;
        G4double sinIncidenceAnglePOut;
    } theIncomingEFunctor(theParticle,theNucleus,theQValueCorrection);

    G4double v = theNucleus->getPotential()->computePotentialEnergy(theParticle);
    if(theParticle->getKineticEnergy()+v-theQValueCorrection<0.) { // Particle entering below 0. Die gracefully
      INCL_DEBUG("Particle " << theParticle->getID() << " is trying to enter below 0" << '\n');
      return false;
    }
    const RootFinder::Solution theSolution = RootFinder::solve(&theIncomingEFunctor, v);
    if(theSolution.success) { // Apply the solution
      theIncomingEFunctor(theSolution.x);
      INCL_DEBUG("Particle successfully entered:\n" << theParticle->print() << '\n');
    } else {
      INCL_WARN("Couldn't compute the potential for incoming particle, root-finding algorithm failed." << '\n');
    }
    return theSolution.success;
  }

}

