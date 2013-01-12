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

/*
 * G4INCLRandom.hh
 *
 *  \date 7 June 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLRANDOM_HH_
#define G4INCLRANDOM_HH_

#include <iostream>
#include <cmath>
#include "G4INCLIRandomGenerator.hh"
#include "G4INCLThreeVector.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLLogger.hh"

namespace G4INCL {

  class Random {
  private:
    Random() {}
    virtual ~Random() {}

  private:
    static IRandomGenerator *theGenerator;

  public:
    /**
     * Set the random number generator implementation to be used globally by INCL.
     *
     * @see G4INCL::IRandomGenerator
     */
    static void setGenerator(G4INCL::IRandomGenerator *aGenerator) {
      if(isInitialized()) {
        ERROR("INCL random number generator already initialized." << std::endl);
      } else {
        theGenerator = aGenerator;
      }
    };

    /**
     * Set the seeds of the current generator.
     *
     */
    static void setSeeds(const SeedVector &sv) {
      theGenerator->setSeeds(sv);
    };

    /**
     * Get the seeds of the current generator.
     *
     */
    static SeedVector getSeeds() {
      return theGenerator->getSeeds();
    };

    /**
     * Generate flat distribution of random numbers.
     */
    static G4double shoot() {return theGenerator->flat(); };

    /**
     * Return a random number in the ]0,1] interval
     */
    static G4double shoot0() {
      G4double r;
      while( (r=shoot()) <= 0. )
        ;
      return r;
    }

    /**
     * Return a random number in the [0,1[ interval
     */
    static G4double shoot1() {
      G4double r;
      while( (r=shoot()) >= 1. )
        ;
      return r;
    }

    /**
     * Generate random numbers using gaussian distribution.
     */
    static G4double gauss(G4double sigma=1.);

    /**
     * Generate isotropically-distributed ThreeVectors of given norm.
     */
    static ThreeVector normVector(G4double norm=1.);

    /**
     * Generate ThreeVectors that are uniformly distributed in a sphere of
     * radius rmax.
     */
    static ThreeVector sphereVector(G4double rmax=1.) {
      return normVector( rmax*Math::pow13(shoot0()) );
    }

    /** \brief Generate Gaussianly-distributed ThreeVectors
     *
     * Generate ThreeVectors that are distributed as a three-dimensional
     * Gaussian of the given sigma.
     */
    static ThreeVector gaussVector(G4double sigma=1.) {
      const G4double sigmax = sigma * Math::oneOverSqrtThree;
      return ThreeVector(gauss(sigmax), gauss(sigmax), gauss(sigmax));
    }

    /**
     * Delete the generator
     */
    static void deleteGenerator() {
      delete theGenerator;
      theGenerator = 0;
    }

    /**
     * Check if the generator is initialized.
     */
    static G4bool isInitialized() {
      if(theGenerator == 0) return false;
      return true;
    };
  };

}

#endif /* G4INCLRANDOM_HH_ */
