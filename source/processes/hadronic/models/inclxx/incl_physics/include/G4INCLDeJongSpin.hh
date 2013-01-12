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

/** \file G4INCLDeJongSpin.hh
 * \brief Simple class implementing De Jong's spin model for nucleus-nucleus
 *        collisions
 *
 * Reference: De Jong, Ignatyuk and Schmidt, Nucl. Phys. A613 (1997) 435-444.
 *
 * \date 2 April 2012
 * \author Davide Mancusi
 */

#ifndef G4INCLDEJONGSPIN_HH_
#define G4INCLDEJONGSPIN_HH_

#include "G4INCLGlobals.hh"
#include "G4INCLRandom.hh"

namespace G4INCL {
  class DeJongSpin {
    public:
      static ThreeVector shoot(const G4int Ap, const G4int Af) {
        return Random::gaussVector(getSpinCutoffParameter(Ap, Af));
      }

    private:
      static G4double getSpinCutoffParameter(const G4int Ap, const G4int Af) {
        const G4double jz2 = jzFactor * Math::pow23((G4double) Ap); // No deformation assumed
        const G4double sigma = jz2 * Af*(Ap-Af)/((G4double)(Ap-1));
        return std::sqrt(sigma);
      }

      const static G4double jzFactor;

      DeJongSpin() {}
      ~DeJongSpin() {}

  };
}

#endif // G4INCLDEJONGSPIN_HH_
