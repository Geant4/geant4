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

/** \file G4INCLIntersection.cc
 * \brief Simple class for computing intersections between a straight line and a sphere.
 *
 * \date 12 December 2011
 * \author Davide Mancusi
 */

#include "G4INCLIntersection.hh"

namespace G4INCL {

  Intersection IntersectionFactory::getTrajectoryIntersection(const ThreeVector &x0, const ThreeVector &v, const G4double r, const G4bool earliest) {
    const G4double scalarVelocity = v.mag();
    ThreeVector velocityUnitVector = v / scalarVelocity;

    ThreeVector positionTransverse = x0 - velocityUnitVector * x0.dot(velocityUnitVector);
    const G4double impactParameter = positionTransverse.mag();

    const G4double r2 = r*r;
    G4double distanceZ2 = r2 - impactParameter * impactParameter;
    if(distanceZ2 < 0.0)
      return Intersection(false, 0.0, ThreeVector());

    const G4double distanceZ = std::sqrt(distanceZ2);
    const ThreeVector position = positionTransverse + velocityUnitVector * (earliest ? -distanceZ : distanceZ);
    const G4double time = (position-x0).dot(velocityUnitVector)/scalarVelocity;
    return Intersection(true, time, position);
  }

}

