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
 * G4INCLRandom.hh
 *
 *  \date 7 June 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLRANDOM_HH_
#define G4INCLRANDOM_HH_

#include <iostream>
#include <cmath>
#include <utility>
#include "G4INCLIRandomGenerator.hh"
#include "G4INCLThreeVector.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {

  namespace Random {
    /**
     * Set the random number generator implementation to be used globally by INCL.
     *
     * @see G4INCL::IRandomGenerator
     */
    void setGenerator(G4INCL::IRandomGenerator *aGenerator);

    /**
     * Set the seeds of the current generator.
     *
     */
    void setSeeds(const SeedVector &sv);

    /**
     * Get the seeds of the current generator.
     *
     */
    SeedVector getSeeds();

    /**
     * Generate flat distribution of random numbers.
     */
    G4double shoot();

    /**
     * Return a random number in the ]0,1] interval
     */
    G4double shoot0();

    /**
     * Return a random number in the [0,1[ interval
     */
    G4double shoot1();

    /**
     * Return a random integer in the [0,n[ interval
     */
    template<typename T> T shootInteger(T n){
      return static_cast<T>(shoot1() * n);
    }

    /**
     * Generate random numbers using gaussian distribution.
     */
    G4double gauss(G4double sigma=1.);

    /**
     * Generate random numbers using gaussian distribution to be used only
     * for correlated pairs
     */
    G4double gaussWithMemory(G4double sigma=1.);

    /**
     * Generate isotropically-distributed ThreeVectors of given norm.
     */
    ThreeVector normVector(G4double norm=1.);

    /**
     * Generate ThreeVectors that are uniformly distributed in a sphere of
     * radius rmax.
     */
    ThreeVector sphereVector(G4double rmax=1.);

    /** \brief Generate Gaussianly-distributed ThreeVectors
     *
     * Generate ThreeVectors that are distributed as a three-dimensional
     * Gaussian of the given sigma.
     */
    ThreeVector gaussVector(G4double sigma=1.);

    /// \brief Generate pairs of correlated Gaussian random numbers
    std::pair<G4double,G4double> correlatedGaussian(const G4double corrCoeff, const G4double x0=0., const G4double sigma=1.);

    /// \brief Generate pairs of correlated uniform random numbers
    std::pair<G4double,G4double> correlatedUniform(const G4double corrCoeff);

    /**
     * Delete the generator
     */
    void deleteGenerator();

    /**
     * Check if the generator is initialized.
     */
    G4bool isInitialized();

#ifdef INCL_COUNT_RND_CALLS
    /// \brief Return the number of calls to the RNG
    unsigned long long getNumberOfCalls();
#endif

    /// \brief Save the status of the random-number generator
    void saveSeeds();

    /// \brief Get the saved status of the random-number generator
    SeedVector getSavedSeeds();

    /// \brief Initialize generator according to a Config object
#ifdef INCLXX_IN_GEANT4_MODE
    void initialize(Config const * const);
#else
    void initialize(Config const * const theConfig);
#endif

    class Adapter {
      public:
        G4int operator()(const G4int n) const;
    };

    Adapter const &getAdapter();
  }

}

#endif /* G4INCLRANDOM_HH_ */
