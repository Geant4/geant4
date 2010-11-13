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
// $Id: G4InclRandomNumbers.hh,v 1.8 2010-11-13 00:08:36 kaitanie Exp $
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

#include "globals.hh"
#include "Randomize.hh"

#ifndef G4InclRandomNumbers_hh
#define G4InclRandomNumbers_hh 1

/**
 * Interface for the random number services.
 */
class G4InclRandomInterface {

public:
  G4InclRandomInterface() {
    this->seed = 1337; // Default seed, this is never actually used.
  }
  G4InclRandomInterface(G4long seed) {
    this->seed = seed;
  }

  virtual ~G4InclRandomInterface() { }

  /**
   * Provide evenly distributed random numbers.
   */
  virtual G4double getRandom() = 0;
  virtual void printSeeds() = 0;
private:
  G4long seed;
};

/**
 * Provides dummy random number generator that always produces one
 * number.
 */
class G4InclDummyRandom : public G4InclRandomInterface {

public:
  G4InclDummyRandom() { }

  ~G4InclDummyRandom() { }

  /**
   * Always return one "random" number. This interface is intended as
   * a debugging tool in order to produce some specific event.
   */
  G4double getRandom() {
    return 0.5;
  }

  void printSeeds() {};
};

/**
 * Interface to the Geant4 random number services.
 */
class G4InclGeant4Random : public G4InclRandomInterface {

public:
  G4InclGeant4Random() { }

  ~G4InclGeant4Random() { }

  /**
   * Always return one "random" number. This interface is intended as
   * a debugging tool in order to produce some specific event.
   */
  G4double getRandom() {
    return G4UniformRand();
  }

  void printSeeds() {
    G4cout <<"Using Geant4 random number generator." << G4endl;
  };
};

#endif
