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
//
// $Id:$
//
// Class G4MTRandGeneral
// 
// Modified version with MT extensions
// of corresponding CLHEP class RandGeneral

// --------------------------------------------------------------------
#ifndef G4MTRandGeneral_h
#define G4MTRandGeneral_h 1

#include <vector>

#include "G4MTHepRandom.hh"

class G4MTRandGeneral : public G4MTHepRandom
{

public:

  G4MTRandGeneral ( const G4double* aProbFunc, 
                    G4int theProbSize, 
                    G4int IntType=0 );
  G4MTRandGeneral ( CLHEP::HepRandomEngine& anEngine,
                    const G4double* aProbFunc, 
                    G4int theProbSize, 
                    G4int IntType=0 );
  G4MTRandGeneral ( CLHEP::HepRandomEngine* anEngine, 
                    const G4double* aProbFunc, 
                    G4int theProbSize, 
                    G4int IntType=0 );
  // These constructors should be used to instantiate a G4MTRandGeneral
  // distribution object defining a local engine for it.
  // The static generator will be skiped by using the non-static methods
  // defined below. In case no engine is specified in the constructor, the
  // default engine used by the static generator is applied.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the G4MTRandGeneral destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the G4MTRandGeneral destructor.
  // The probability distribution function (Pdf) must be provided by the user
  // as an array of positive real number. The array size must also be
  // provided. The Pdf doesn't need to be normalized to 1. 
  // if IntType = 0 ( default value ) a uniform random number is
  // generated using the engine. The uniform number is then transformed
  // to the user's distribution using the cumulative probability
  // distribution constructed from his histogram. The cumulative
  // distribution is inverted using a binary search for the nearest
  // bin boundary and a linear interpolation within the
  // bin. G4MTRandGeneral therefore generates a constant density within
  // each bin.
  // if IntType = 1 no interpolation is performed and the result is a
  // discrete distribution.

  virtual ~G4MTRandGeneral();
  // Destructor

  // Methods to shoot random values using the static generator
  // N.B.: The methods are NOT static since they use nonstatic members
  // theIntegralPdf & nBins

  inline G4double shoot();

  inline void shootArray ( const G4int size, G4double* vect);

  //  Methods to shoot random values using a given engine
  //  by-passing the static generator.

  G4double shoot( CLHEP::HepRandomEngine* anEngine );

  void shootArray ( CLHEP::HepRandomEngine* anEngine, const G4int size,
                    G4double* vect );
    
  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  G4double fire();

  void fireArray ( const G4int size, G4double* vect);

  G4double operator()();

private:

  CLHEP::HepRandomEngine* localEngine;
  G4bool deleteEngine;
  std::vector<G4double> theIntegralPdf;
  G4int nBins;
  G4double oneOverNbins;
  G4int InterpolationType;

  // Private methods to factor out replicated implementation sections
  void prepareTable(const G4double* aProbFunc);
  void useFlatDistribution();
  G4double mapRandom(G4double rand) const;

};

#include "G4MTRandGeneral.icc"

#endif
