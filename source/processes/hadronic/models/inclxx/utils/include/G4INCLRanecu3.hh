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

/*
 * \file G4INCLRanecu3.hh
 *
 *  \date 25 October 2013
 * \author Davide Mancusi
 */

#ifndef G4INCLRanecu3_HH_
#define G4INCLRanecu3_HH_

#include "G4INCLIRandomGenerator.hh"
// #include <cassert>

namespace G4INCL {

  /** \brief Extended Ranecu-type RNG class
   *
   * This generator implements the C++ version of the generator suggested by
   * Badal and Sempau, Comp. Phys. Comm. 175 (2006) 440. It uses three 32-bit
   * seeds and has a period of ~5E27.
   */
  class Ranecu3 : public G4INCL::IRandomGenerator {
  public:
    Ranecu3();
    Ranecu3(const Random::SeedVector &sv);
    virtual ~Ranecu3();

    Random::SeedVector getSeeds() {
      Random::SeedVector sv;
      sv.push_back(iseed1);
      sv.push_back(iseed2);
      sv.push_back(iseed3);
      return sv;
    }

    void setSeeds(const Random::SeedVector &sv) {
// assert(sv.size()>=3);
      iseed1 = sv[0];
      iseed2 = sv[1];
      iseed3 = sv[2];
    }

    G4double flat();

  private:
    G4int iseed1, iseed2, iseed3;
    G4int i1, i2, i3, iz;
    const G4double uscale;
    const G4int m1, m2, m3;
    const G4int a1, a2, a3;
    const G4int q1, q2, q3;
    const G4int r1, r2, r3;
  };

}

#endif /* G4INCLRanecu3_HH_ */
