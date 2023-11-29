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

/** \file G4INCLCoulombNone.cc
 * \brief Placeholder class for no Coulomb distortion.
 *
 * \date 14 February 2011
 * \author Davide Mancusi
 */

#include "G4INCLCoulombNone.hh"
#include "G4INCLIntersection.hh"

namespace G4INCL {

  ParticleEntryAvatar *CoulombNone::bringToSurface(Particle * const p, Nucleus * const n) const {
    Intersection intersection = IntersectionFactory::getEarlierTrajectoryIntersection(p->getPosition(), p->getPropagationVelocity(), n->getUniverseRadius());
    if(intersection.exists) { // If the particle enters the nucleus
      p->setPosition(intersection.position);
      return new ParticleEntryAvatar(0.0, n, p);
    } else // If the particle does NOT enter the nucleus
      return NULL;
  }

  IAvatarList CoulombNone::bringToSurface(Cluster * const c, Nucleus * const n) const {
    // The avatar list that we will return
    IAvatarList theAvatarList;

    // Loop over the particles in the cluster
    ParticleList const &projectiles = c->getParticles();
    std::list<Intersection> theIntersections;
    G4double theFirstEntryTime = 1E+60; // a large time
    G4int theFirstID = 0;
    for(ParticleIter p=projectiles.begin(), e=projectiles.end(); p!=e; ++p) {
      // Check if the particle enters the nucleus
      Intersection intersection(IntersectionFactory::getEarlierTrajectoryIntersection(
            (*p)->getPosition(),
            (*p)->getPropagationVelocity(),
            n->getUniverseRadius()));
      // Store the intersections
      theIntersections.push_back(intersection);
      if(intersection.exists) {
        // Position the particle at the entry point
        (*p)->setPosition(intersection.position);

        // Keep track of the first entering particle
        if(intersection.time < theFirstEntryTime) {
          theFirstEntryTime = intersection.time;
          theFirstID = (G4int)(*p)->getID();
        }
      }
    }

    std::list<Intersection>::const_iterator intIter = theIntersections.begin();
    for(ParticleIter p=projectiles.begin(), e=projectiles.end(); p!=e; ++p, ++intIter) {

      if((*intIter).exists) {
        // If the particle enters the nucleus, generate a ParticleEntryAvatar
        // for it and add it to the list of avatars that we will return
        if((*p)->getID() == theFirstID) {
          // The first particle always enters exactly at t=0 (in order to
          // avoid negative entry times due to rounding)
          theAvatarList.push_back(new ParticleEntryAvatar(0.0, n, *p));
        } else
          theAvatarList.push_back(new ParticleEntryAvatar(intIter->time - theFirstEntryTime, n, *p));
      }

   }

    return theAvatarList;
  }

}
