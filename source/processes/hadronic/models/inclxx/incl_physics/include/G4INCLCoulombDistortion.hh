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

/** \file G4INCLCoulombDistortion.hh
 * \brief Static class for selecting Coulomb distortion.
 *
 * \date 14 February 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLCOULOMBDISTORTION_HH_
#define G4INCLCOULOMBDISTORTION_HH_

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLICoulomb.hh"

namespace G4INCL {

  /**
   * Coulomb distortion
   */
  class CoulombDistortion {
  public:
    /** \brief Modify the momentum of an incoming particle and position it on
     *         the surface of the nucleus.
     *
     * This method places Particle p on the surface of Nucleus n and modifies
     * the direction of its momentum to be tangent to the Coulomb trajectory in
     * that point.
     *
     * The input particle has to be prepared with its asymptotic momentum. Its
     * position is used only for the purpose of computing the asymptotic impact
     * parameter; in other words, this method only uses the components of the
     * position that are perpendicular to the momentum. The remaining component
     * is not used, and can be set to any value.
     *
     * This method returns a ParticleEntry avatar for the projectile.
     *
     * \param p incoming particle
     * \param n distorting nucleus
     * \return the ParticleEntryAvatar for the projectile particle
     **/
    static ParticleEntryAvatar *bringToSurface(Particle *p, Nucleus * const n) {
      return theCoulomb->bringToSurface(p, n);
    }

    /** \brief Modify the momentum of an incoming cluster and position it on
     *         the surface of the target.
     *
     * Same as the Particle-based bringToSurface method, but for incoming heavy
     * ions.
     *
     * This method returns a list of ParticleEntry avatars for the participant
     * nucleons
     *
     * \param c incoming heavy ion
     * \param n distorting nucleus
     * \return a list of ParticleEntryAvatars
     **/
    static IAvatarList bringToSurface(Cluster * const c, Nucleus * const n) {
      return theCoulomb->bringToSurface(c, n);
    }

    /** \brief Modify the momentum of an outgoing particle. */
    static void distortOut(ParticleList const &pL, Nucleus const * const n) {
      theCoulomb->distortOut(pL, n);
    }

    /** \brief Return the maximum impact parameter for Coulomb-distorted
     *         trajectories. **/
    static G4double maxImpactParameter(ParticleSpecies const &p, const G4double kinE, Nucleus const * const n) {
      return theCoulomb->maxImpactParameter(p, kinE, n);
    }

    /** \brief Return the maximum impact parameter for Coulomb-distorted
     *         trajectories. **/
    static G4double maxImpactParameter(Particle const * const p, Nucleus const * const n) {
      return maxImpactParameter(p->getSpecies(), p->getKineticEnergy(), n);
    }

    /** \brief Set the Coulomb-distortion algorithm. */
    static void setCoulomb(ICoulomb * const coulomb) { theCoulomb = coulomb; }

    /** \brief Delete the Coulomb-distortion object. */
    static void deleteCoulomb() {
      delete theCoulomb;
      theCoulomb = 0;
    }

  protected:
    CoulombDistortion() {}
    ~CoulombDistortion() {}

  private:
    static ICoulomb *theCoulomb;

  };
}

#endif /* G4INCLCOULOMBDISTORTION_HH_ */
