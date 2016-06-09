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

/** \file G4INCLCoulombNonRelativistic.hh
 * \brief Class for non-relativistic Coulomb distortion.
 *
 * Created on: 14 February 2011
 *     Author: Davide Mancusi
 */

#ifndef G4INCLCOULOMBNONRELATIVISTIC_HH_
#define G4INCLCOULOMBNONRELATIVISTIC_HH_

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLICoulomb.hh"

namespace G4INCL {

  class CoulombNonRelativistic : public ICoulomb {
  public:
    CoulombNonRelativistic() {}
    virtual ~CoulombNonRelativistic() {}

    /** \brief Modify the momentum of the particle and position it on the
     *         surface of the nucleus.
     *
     * This method performs non-relativistic distortion.
     *
     * \param p incoming particle
     * \param n distorting nucleus
     **/
    void bringToSurface(Particle * const p, Nucleus const * const n) const;

    /** \brief Modify the momenta of the outgoing particles.
     *
     * This method performs non-relativistic distortion.
     *
     * \param pL list of outgoing particles
     * \param n distorting nucleus
     */
    void distortOut(ParticleList const &pL, Nucleus const * const n) const;

    /** \brief Return the maximum impact parameter for Coulomb-distorted
     *         trajectories. **/
    G4double maxImpactParameter(Particle const * const p, Nucleus const *
        const n) const;

  private:
    /** \brief Return the "Coulomb factor". */
    G4double coulombFactor(Particle const * const p, Nucleus const * const n) const {
      return eSquared * p->getZ() * n->getZ() / p->getKineticEnergy();
    }
  };
}

#endif /* G4INCLCOULOMBNONRELATIVISTIC_HH_ */
