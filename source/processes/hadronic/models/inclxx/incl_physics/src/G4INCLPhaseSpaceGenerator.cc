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

#include "G4INCLPhaseSpaceGenerator.hh"
#include "G4INCLPhaseSpaceKopylov.hh"
#include "G4INCLPhaseSpaceRauboldLynch.hh"

namespace G4INCL {

  namespace {
    G4ThreadLocal IPhaseSpaceGenerator *thePhaseSpaceGenerator;

    G4ThreadLocal Particle *biasMe;

    /** \brief Actually perform the biasing
     *
     * \param particles list of particles to bias
     * \param pInVec momentum of the particle to be biased before the collision
     * \param slope the parameter \f$B\f$ in \f$\exp(B\cdot t)\f$
     */
    void bias(ParticleList &particles, const ThreeVector &pInVec, const G4double slope) {
      const G4double pIn = pInVec.mag();
      const ThreeVector collisionAxis = pInVec/pIn;
      const ThreeVector pMomVec = biasMe->getMomentum();
      const G4double pMom = pMomVec.mag();
      const G4double pMomCosAng = pMomVec.dot(collisionAxis)/pMom;
      const G4double pMomAng = Math::arcCos(pMomCosAng); // Angle between the original axis of the dominant particle and is new one after generate

      // compute the target angle for the biasing
      // it is drawn from a exp(Bt) distribution
      const G4double cosAngSlope = 2e-6 * slope * pIn * pMom;
      const G4double cosAng = 1. + std::log(1. - Random::shoot()*(1.-std::exp(-2.*cosAngSlope)))/cosAngSlope;
      const G4double ang = Math::arcCos(cosAng);

      // compute the rotation angle
      const G4double rotationAngle = ang - pMomAng;

      // generate the rotation axis; it is perpendicular to collisionAxis and
      // pMomVec
      ThreeVector rotationAxis;
      if(pMomAng>1E-10) {
        rotationAxis = collisionAxis.vector(pMomVec);
        const G4double axisLength = rotationAxis.mag();
        const G4double oneOverLength = 1./axisLength;
        rotationAxis *= oneOverLength;
      } else {
        // need to jump through some hoops if collisionAxis is nearly aligned
        // with pMomVec
        rotationAxis = collisionAxis.anyOrthogonal();
      }

      // apply the rotation
      particles.rotateMomentum(rotationAngle, rotationAxis);
    }

  }

  namespace PhaseSpaceGenerator {
    void generate(const G4double sqrtS, ParticleList &particles) {
      return thePhaseSpaceGenerator->generate(sqrtS, particles);
    }

    void generateBiased(const G4double sqrtS, ParticleList &particles, const size_t index, const G4double slope) {
// assert(index<particles.size());
      // store the incoming momentum of particle[index]; it will be used to
      // compute t when biasing
      biasMe = particles[index];
      const ThreeVector pInVec = biasMe->getMomentum();
      generate(sqrtS, particles);
      bias(particles, pInVec, slope);
    }

    void setPhaseSpaceGenerator(IPhaseSpaceGenerator *g) {
      thePhaseSpaceGenerator = g;
    }

    IPhaseSpaceGenerator *getPhaseSpaceGenerator() {
      return thePhaseSpaceGenerator;
    }

    void deletePhaseSpaceGenerator() {
      delete thePhaseSpaceGenerator;
      thePhaseSpaceGenerator = NULL;
    }

    void initialize(Config const * const theConfig) {
      PhaseSpaceGeneratorType psg = theConfig->getPhaseSpaceGeneratorType();
      if(psg==RauboldLynchType)
        setPhaseSpaceGenerator(new PhaseSpaceRauboldLynch);
      else if(psg==KopylovType)
        setPhaseSpaceGenerator(new PhaseSpaceKopylov);
      else
        setPhaseSpaceGenerator(NULL);
    }
  }
}
