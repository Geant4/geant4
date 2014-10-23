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

/** \file G4INCLCoulombDistortion.cc
 * \brief Static class for selecting Coulomb distortion.
 *
 * \date 14 February 2011
 * \author Davide Mancusi
 */

#include "G4INCLCoulombDistortion.hh"
#include "G4INCLCoulombNone.hh"
#include "G4INCLCoulombNonRelativistic.hh"

namespace G4INCL {

  namespace CoulombDistortion {

    namespace {
      G4ThreadLocal ICoulomb* theCoulomb = 0;
    }

    ParticleEntryAvatar *bringToSurface(Particle *p, Nucleus * const n) {
      return theCoulomb->bringToSurface(p, n);
    }

    IAvatarList bringToSurface(Cluster * const c, Nucleus * const n) {
      return theCoulomb->bringToSurface(c, n);
    }

    void distortOut(ParticleList const &pL, Nucleus const * const n) {
      theCoulomb->distortOut(pL, n);
    }

    G4double maxImpactParameter(ParticleSpecies const &p, const G4double kinE, Nucleus const * const n) {
      return theCoulomb->maxImpactParameter(p, kinE, n);
    }

    G4double maxImpactParameter(Particle const * const p, Nucleus const * const n) {
      return maxImpactParameter(p->getSpecies(), p->getKineticEnergy(), n);
    }

    void setCoulomb(ICoulomb * const coulomb) { theCoulomb = coulomb; }

    void deleteCoulomb() {
      delete theCoulomb;
      theCoulomb = 0;
    }

    void initialize(Config const * const theConfig) {
      CoulombType coulombType = theConfig->getCoulombType();
      if(coulombType == NonRelativisticCoulomb)
        setCoulomb(new CoulombNonRelativistic);
      else if(coulombType == NoCoulomb)
        setCoulomb(new CoulombNone);
      else
        setCoulomb(NULL);
    }
  }

}
