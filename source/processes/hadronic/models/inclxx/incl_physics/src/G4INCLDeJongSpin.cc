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

/** \file G4INCLDeJongSpin.cc
 * \brief Simple class implementing De Jong's spin model for nucleus-nucleus
 *        collisions
 *
 * Reference: De Jong, Ignatyuk and Schmidt, Nucl. Phys. A613 (1997) 435-444.
 *
 * \date 2 April 2012
 * \author Davide Mancusi
 */

#include "G4INCLGlobals.hh"
#include "G4INCLRandom.hh"
#include "G4INCLDeJongSpin.hh"

namespace G4INCL {

  namespace DeJongSpin {

    namespace {

      const G4double jzFactor = PhysicalConstants::hcSquared * 0.16;

      G4double getSpinCutoffParameter(const G4int Ap, const G4int Af) {
        const G4double jz2 = jzFactor * Math::pow23((G4double) Ap); // No deformation assumed
        const G4double sigma = jz2 * Af*(Ap-Af)/((G4double)(Ap-1));
        return std::sqrt(sigma);
      }

    }

    ThreeVector shoot(const G4int Ap, const G4int Af) {
      return Random::gaussVector(getSpinCutoffParameter(Ap, Af));
    }
  }
}
