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
#include "G4INCLConfig.hh"

namespace G4INCL {

  /**
   * Coulomb distortion
   */
  namespace CoulombDistortion {
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
    ParticleEntryAvatar *bringToSurface(Particle *p, Nucleus * const n);

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
    IAvatarList bringToSurface(Cluster * const c, Nucleus * const n);

    /** \brief Modify the momentum of an outgoing particle. */
    void distortOut(ParticleList const &pL, Nucleus const * const n);

    /** \brief Return the maximum impact parameter for Coulomb-distorted
     *         trajectories. **/
    G4double maxImpactParameter(ParticleSpecies const &p, const G4double kinE, Nucleus const * const n);

    /** \brief Return the maximum impact parameter for Coulomb-distorted
     *         trajectories. **/
    G4double maxImpactParameter(Particle const * const p, Nucleus const * const n);

    /** \brief Set the Coulomb-distortion algorithm. */
    void setCoulomb(ICoulomb * const coulomb);

    /** \brief Delete the Coulomb-distortion object. */
    void deleteCoulomb();

    /// \brief Initialize the Coulomb-distortion algorithm
    void initialize(Config const * const theConfig);

  }
}

#endif /* G4INCLCOULOMBDISTORTION_HH_ */
