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

#ifndef G4INCLGeant4Random_hh
#define G4INCLGeant4Random_hh 1

#ifdef INCLXX_IN_GEANT4_MODE

#include "globals.hh"
#include "Randomize.hh"

#include "G4INCLLogger.hh"
#include "G4INCLIRandomGenerator.hh"

namespace G4INCL {

  class Geant4RandomGenerator : public G4INCL::IRandomGenerator {
  public:
    Geant4RandomGenerator() {};
    Geant4RandomGenerator(const Random::SeedVector &) {};
    virtual ~Geant4RandomGenerator() {};

    Random::SeedVector getSeeds() {
      INCL_WARN("getSeeds not supported.");
      Random::SeedVector sv;
      return sv;
    }

    void setSeeds(const Random::SeedVector &) {
      INCL_WARN("setSeeds not supported.");
    }

    G4double flat() {
      return G4UniformRand();
    }
  };

}

#endif // INCLXX_IN_GEANT4_MODE

#endif
