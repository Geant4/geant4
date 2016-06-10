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

/** \file G4INCLIntersection.hh
 * \brief Simple class for computing intersections between a straight line and a sphere.
 *
 * \date 12 December 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLINTERSECTION_HH
#define G4INCLINTERSECTION_HH 1

#include "G4INCLThreeVector.hh"
#include <utility>

namespace G4INCL {

  /** \brief Intersection-point structure
   *
   * The structure contains the time and position of the intersection point
   * of a trajectory with a surface, if it exists.
   */
  struct Intersection {
    Intersection(const G4bool e, const G4double t, const ThreeVector &p) :
      exists(e), time(t), position(p) {}
    G4bool exists;
    G4double time;
    ThreeVector position;
  };

  namespace IntersectionFactory {

    /** \brief Compute the first intersection of a straight particle
     *         trajectory with a sphere.
     *
     * \param x0 the starting position of the trajectory
     * \param p the trajectory direction
     * \param r the radius of the sphere (centred in the origin)
     * \return an Intersection. The G4bool is true if an intersection exists,
     *         in which case its position is stored in the ThreeVector and
     *         its time in the G4double.
     */
    Intersection getEarlierTrajectoryIntersection(const ThreeVector &x0, const ThreeVector &p, const G4double r);
    /** \brief Compute the second intersection of a straight particle
     *         trajectory with a sphere.
     *
     * \param x0 the starting position of the trajectory
     * \param p the trajectory direction
     * \param r the radius of the sphere (centred in the origin)
     * \return an Intersection. The G4bool is true if an intersection exists,
     *         in which case its position is stored in the ThreeVector and
     *         its time in the G4double.
     */
    Intersection getLaterTrajectoryIntersection(const ThreeVector &x0, const ThreeVector &p, const G4double r);
    /** \brief Compute both intersections of a straight particle
     *         trajectory with a sphere.
     *
     * \param x0 the starting position of the trajectory
     * \param p the trajectory direction
     * \param r the radius of the sphere (centred in the origin)
     * \return an Intersection. The G4bool is true if an intersection exists,
     *         in which case its position is stored in the ThreeVector and
     *         its time in the G4double.
     */
    std::pair<Intersection,Intersection> getTrajectoryIntersections(const ThreeVector &x0, const ThreeVector &p, const G4double r);

    namespace {

      Intersection getTrajectoryIntersection(const ThreeVector &x0, const ThreeVector &v, const G4double r, const G4bool earliest) {
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

    inline Intersection getEarlierTrajectoryIntersection(const ThreeVector &x0, const ThreeVector &p, const G4double r) {
      return getTrajectoryIntersection(x0, p, r, true);
    }

    inline Intersection getLaterTrajectoryIntersection(const ThreeVector &x0, const ThreeVector &p, const G4double r) {
      return getTrajectoryIntersection(x0, p, r, false);
    }

    inline std::pair<Intersection,Intersection> getTrajectoryIntersections(const ThreeVector &x0, const ThreeVector &p, const G4double r) {
      return std::make_pair(
                            getTrajectoryIntersection(x0, p, r, true),
                            getTrajectoryIntersection(x0, p, r, false)
                           );
    }
  }

}

#endif /* G4INCLINTERSECTION_HH */
