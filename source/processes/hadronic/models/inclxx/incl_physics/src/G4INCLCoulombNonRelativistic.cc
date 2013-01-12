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

/** \file G4INCLCoulombNonRelativistic.cc
 * \brief Class for non-relativistic Coulomb distortion.
 *
 * \date 14 February 2011
 * \author Davide Mancusi
 */

#include "G4INCLCoulombNonRelativistic.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  ParticleEntryAvatar *CoulombNonRelativistic::bringToSurface(Particle * const p, Nucleus * const n) const {
    // No distortion for neutral particles
    if(p->getZ()!=0) {
      const G4bool success = coulombDeviation(p, n);
      if(!success) // transparent
        return NULL;
    }

    // Rely on the CoulombNone slave to compute the straight-line intersection
    // and actually bring the particle to the surface of the nucleus
    return theCoulombNoneSlave.bringToSurface(p,n);
  }

  IAvatarList CoulombNonRelativistic::bringToSurface(Cluster * const c, Nucleus * const n) const {
    // Neutral clusters?!
// assert(c->getZ()>0);

    // Perform the actual Coulomb deviation
    const G4bool success = coulombDeviation(c, n);
    if(!success) {
      return IAvatarList();
    }

    // Rely on the CoulombNone slave to compute the straight-line intersection
    // and actually bring the particle to the surface of the nucleus
    return theCoulombNoneSlave.bringToSurface(c,n);
  }

  void CoulombNonRelativistic::distortOut(ParticleList const &pL,
      Nucleus const * const nucleus) const {

    for(ParticleIter particle=pL.begin(); particle!=pL.end(); ++particle) {

      const G4int Z = (*particle)->getZ();
      if(Z == 0) continue;

      const G4double tcos=1.-0.000001;

      const G4double et1 = PhysicalConstants::eSquared * nucleus->getZ();
      const G4double transmissionRadius =
        nucleus->getDensity()->getTransmissionRadius(*particle);

      const ThreeVector position = (*particle)->getPosition();
      ThreeVector momentum = (*particle)->getMomentum();
      const G4double r = position.mag();
      const G4double p = momentum.mag();
      const G4double cosTheta = position.dot(momentum)/(r*p);
      if(cosTheta < 0.999999) {
        const G4double sinTheta = std::sqrt(1.-cosTheta*cosTheta);
        const G4double eta = et1 * Z / (*particle)->getKineticEnergy();
        if(eta > transmissionRadius-0.0001) {
          // If below the Coulomb barrier, radial emission:
          momentum = position * (p/r);
          (*particle)->setMomentum(momentum);
        } else {
          const G4double b0 = 0.5 * (eta + std::sqrt(eta*eta +
                4. * std::pow(transmissionRadius*sinTheta,2)
                * (1.-eta/transmissionRadius)));
          const G4double bInf = std::sqrt(b0*(b0-eta));
          const G4double thr = std::atan(eta/(2.*bInf));
          G4double uTemp = (1.-b0/transmissionRadius) * std::sin(thr) +
            b0/transmissionRadius;      
          if(uTemp>tcos) uTemp=tcos;
          const G4double thd = std::acos(cosTheta)-Math::piOverTwo + thr +
            std::acos(uTemp);
          const G4double c1 = std::sin(thd)*cosTheta/sinTheta + std::cos(thd);
          const G4double c2 = -p*std::sin(thd)/(r*sinTheta);
          const ThreeVector newMomentum = momentum*c1 + position*c2;
          (*particle)->setMomentum(newMomentum);
        }
      }
    }
  }

  G4double CoulombNonRelativistic::maxImpactParameter(ParticleSpecies const &p, const G4double kinE,
      Nucleus const * const n) const {
    G4double theMaxImpactParameter = maxImpactParameterParticle(p, kinE, n);
    if(theMaxImpactParameter <= 0.)
      return 0.;
    if(p.theType == Composite)
      theMaxImpactParameter +=  2.*ParticleTable::getNuclearRadius(p.theA, p.theZ);
    return theMaxImpactParameter;
  }

  G4double CoulombNonRelativistic::maxImpactParameterParticle(ParticleSpecies const &p, const G4double kinE,
      Nucleus const * const n) const {
    const G4double theMinimumDistance = minimumDistance(p, kinE, n);
    const G4double rMax = n->getCoulombRadius(p);
    const G4double theMaxImpactParameterSquared = rMax*(rMax-theMinimumDistance);
    if(theMaxImpactParameterSquared<=0.)
      return 0.;
    G4double theMaxImpactParameter = std::sqrt(theMaxImpactParameterSquared);
    return theMaxImpactParameter;
  }

  G4bool CoulombNonRelativistic::coulombDeviation(Particle * const p, Nucleus const * const n) const {
    // Determine the rotation angle and the new impact parameter
    ThreeVector positionTransverse = p->getTransversePosition();
    const G4double impactParameter = positionTransverse.mag();

    // Some useful variables
    const G4double theMinimumDistance = minimumDistance(p, n);
    // deltaTheta2 = (pi - Rutherford scattering angle)/2
    const G4double deltaTheta2 = std::atan(2.*impactParameter/theMinimumDistance);
    const G4double eccentricity = 1./std::cos(deltaTheta2);

    G4double newImpactParameter, alpha; // Parameters that must be determined by the deviation

    ParticleSpecies aSpecies = p->getSpecies();
    G4double kineticEnergy = p->getKineticEnergy();
    // Note that in the following call to maxImpactParameter we are not
    // interested in the size of the cluster. This is why we call
    // maxImpactParameterParticle.
    if(impactParameter>maxImpactParameterParticle(aSpecies, kineticEnergy, n)) {
      // This should happen only for composite particles, whose trajectory can
      // geometrically miss the nucleus but still trigger a cascade because of
      // the finite extension of the projectile.
      // In this case, the sphere radius is the minimum distance of approach
      // and the kinematics is very simple.
      newImpactParameter = 0.5 * theMinimumDistance * (1.+eccentricity); // the minimum distance of approach
      alpha = Math::piOverTwo - deltaTheta2; // half the Rutherford scattering angle
    } else {
      // The particle trajectory intersects the Coulomb sphere

      // Compute the entrance angle
      const G4double radius = n->getCoulombRadius(p->getSpecies());
      G4double argument = -(1. + 2.*impactParameter*impactParameter/(radius*theMinimumDistance))
        / eccentricity;
      const G4double thetaIn = Math::twoPi - std::acos(argument) - deltaTheta2;

      // Velocity angle at the entrance point
      alpha = std::atan((1+std::cos(thetaIn))
        / (std::sqrt(eccentricity*eccentricity-1.) - std::sin(thetaIn)));
      // New impact parameter
      newImpactParameter = radius * std::sin(thetaIn - alpha);
    }

    // Modify the impact parameter of the particle
    positionTransverse *= newImpactParameter/positionTransverse.mag();
    const ThreeVector theNewPosition = p->getLongitudinalPosition() + positionTransverse;
    p->setPosition(theNewPosition);

    // Determine the rotation axis for the incoming particle
    const ThreeVector &momentum = p->getMomentum();
    ThreeVector rotationAxis = momentum.vector(positionTransverse);
    const G4double axisLength = rotationAxis.mag();
    // Apply the rotation
    if(axisLength>1E-20) {
      rotationAxis /= axisLength;
      p->rotate(alpha, rotationAxis);
    }

    return true;
  }

}
