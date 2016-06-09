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

/** \file G4INCLCoulombNonRelativistic.cc
 * \brief Class for non-relativistic Coulomb distortion.
 *
 * Created on: 14 February 2011
 *     Author: Davide Mancusi
 */

#include "G4INCLCoulombNonRelativistic.hh"
#include "G4INCLGlobals.hh"

namespace G4INCL {

  void CoulombNonRelativistic::bringToSurface(Particle * const p, Nucleus const * const n) const {
    ThreeVector momentumUnitVector = p->getMomentum();
    momentumUnitVector /= momentumUnitVector.mag();

    ThreeVector positionTransverse = p->getTransversePosition();
    const G4double impactParameter = positionTransverse.mag();

    const G4double radius = n->getSurfaceRadius(p);

    // No distortion for neutral particles
    G4double newImpactParameter;
    G4double alpha;
    if(p->getZ()==0) {
      newImpactParameter = impactParameter;
      alpha = 0.;
    } else {
      const G4double theCoulombFactor = coulombFactor(p, n);
      const G4double thrs2 = std::atan(theCoulombFactor/(2*impactParameter));
      const G4double eccentricity = -1./std::sin(thrs2);
      const G4double bMin = 0.5 * (
          theCoulombFactor +
          std::sqrt(theCoulombFactor*theCoulombFactor +
            4*impactParameter*impactParameter)
          );
      const G4double phyp = (1.+eccentricity) * bMin;
      const G4double thetaMax = std::acos((phyp/radius - 1.)/eccentricity);
      newImpactParameter = radius * std::cos(thrs2 + thetaMax);
      const G4double psi = std::atan( (1.+eccentricity*std::cos(thetaMax)) /
          (eccentricity*std::sin(thetaMax)));
      alpha = psi - Math::piOverTwo + thrs2 + thetaMax;
    }
    const G4double distanceZ2 = radius*radius - newImpactParameter*newImpactParameter;
    const G4double distanceZ = (distanceZ2>0. ? std::sqrt(distanceZ2) : 0.);

    positionTransverse *= newImpactParameter/impactParameter;

    const ThreeVector position = positionTransverse - momentumUnitVector *
      distanceZ;
    p->setPosition(position);

    positionTransverse /= positionTransverse.mag();
    const G4double momentum = p->getMomentum().mag();
    const ThreeVector newMomentum = p->getMomentum() * std::cos(alpha) +
      positionTransverse * (std::sin(alpha) * momentum);

    p->setMomentum(newMomentum);

  }

  void CoulombNonRelativistic::distortOut(ParticleList const &pL,
      Nucleus const * const nucleus) const {

    for(ParticleIter particle=pL.begin(); particle!=pL.end(); ++particle) {

      const G4int Z = (*particle)->getZ();
      if(Z == 0) continue;

      const G4double tcos=1.-0.000001;

      const G4double et1 = eSquared * nucleus->getZ();
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

  G4double CoulombNonRelativistic::maxImpactParameter(Particle const * const p,
      Nucleus const * const n) const {
    const G4double theCoulombFactor = coulombFactor(p, n);
    const G4double rMax = n->getSurfaceRadius(p);
    const G4double theMaxImpactParameterSquared = rMax*(rMax-theCoulombFactor);
    return (theMaxImpactParameterSquared>0. ?
        std::sqrt(theMaxImpactParameterSquared) : 0.);
  }
}
