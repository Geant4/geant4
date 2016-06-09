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

/** \file G4INCLCoulombDistortion.hh
 * \brief Static class for selecting Coulomb distortion.
 *
 * Created on: 14 February 2011
 *     Author: Davide Mancusi
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
    /** \brief Modify the momentum of an incoming particle.
     *
     * This method places Particle p on the surface of Nucleus n and modifies
     * the direction of its momentum to be tangent to the Coulomb trajectory in
     * that poG4int.
     *
     * The input particle has to be prepared with its asymptotic momentum. Its
     * position is used only for the purpose of computing the asymptotic impact
     * parameter; in other words, this method only uses the components of the
     * position that are perpendicular to the momentum. The remaining component
     * is not used, and can be set to any value.
     *
     * \param p incoming particle
     * \param n distorting nucleus
     **/
    static void bringToSurface(Particle *p, Nucleus const * const n) {
      theCoulomb->bringToSurface(p, n);
    }

    /** \brief Modify the momentum of an outgoing particle. */
    static void distortOut(ParticleList const &pL, Nucleus const * const n) {
      theCoulomb->distortOut(pL, n);
    }

    /** \brief Return the maximum impact parameter for Coulomb-distorted
     *         trajectories. **/
    static G4double maxImpactParameter(Particle const * const p, Nucleus const * const n) {
      return theCoulomb->maxImpactParameter(p, n);
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
