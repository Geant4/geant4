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

/** \file G4INCLCoulombNone.hh
 * \brief Placeholder class for no Coulomb distortion.
 *
 * \date 14 February 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLCOULOMBNONE_HH_
#define G4INCLCOULOMBNONE_HH_

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLICoulomb.hh"
#include <utility>

namespace G4INCL {

  class CoulombNone : public ICoulomb {

    public:
    CoulombNone() {}
    virtual ~CoulombNone() {}

    /** \brief Position the particle on the surface of the nucleus.
     *
     * This method does not perform any distortion.
     *
     * \param p incoming particle
     * \param n distorting nucleus
     **/
    ParticleEntryAvatar *bringToSurface(Particle * const p, Nucleus * const n) const;

    /** \brief Position the cluster on the surface of the nucleus.
     *
     * This method does not perform any distortion.
     *
     * \param c incoming cluster
     * \param n distorting nucleus
     **/
    IAvatarList bringToSurface(Cluster * const c, Nucleus * const n) const;

    /** \brief Modify the momenta of the outgoing particles.
     *
     * This method does not perform any distortion.
     */
    void distortOut(ParticleList const & /* pL */, Nucleus const * const /* n */) const {}

    /** \brief Return the maximum impact parameter for Coulomb-distorted
     *         trajectories. **/
    G4double maxImpactParameter(ParticleSpecies const &p, const G4double /*kinE*/, Nucleus const *
        const n) const {
      if(p.theType == Composite)
        return 2.*ParticleTable::getLargestNuclearRadius(p.theA, p.theZ)
          + n->getUniverseRadius();
      else
        return n->getUniverseRadius();
    }

  };
}

#endif /* G4INCLCOULOMBNONE_HH_ */
