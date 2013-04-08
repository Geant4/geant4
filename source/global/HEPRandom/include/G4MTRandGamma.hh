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
// Class G4MTRandGamma
// 
// Modified version with MT extensions
// of corresponding CLHEP class RandGamma

// --------------------------------------------------------------------
#ifndef G4MTRandGamma_h
#define G4MTRandGamma_h 1

#include "G4MTHepRandom.hh"

class G4MTRandGamma : public G4MTHepRandom
{

public:

  inline G4MTRandGamma ( CLHEP::HepRandomEngine& anEngine, G4double k=1.0,
                                                G4double lambda=1.0 );
  inline G4MTRandGamma ( CLHEP::HepRandomEngine* anEngine, G4double k=1.0, 
                                                G4double lambda=1.0 );
  // These constructors should be used to instantiate a G4MTRandGamma
  // distribution object defining a local engine for it.
  // The static generator will be skipped using the non-static methods
  // defined below.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the G4MTRandGamma destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the G4MTRandGamma destructor.

  virtual ~G4MTRandGamma();
  // Destructor

  // Static methods to shoot random values using the static generator

  static inline G4double shoot();

  static G4double shoot( G4double k, G4double lambda );

  static void shootArray ( const G4int size, G4double* vect,
                            G4double k=1.0, G4double lambda=1.0 );

  //  Static methods to shoot random values using a given engine
  //  by-passing the static generator.

  static inline G4double shoot( CLHEP::HepRandomEngine* anEngine );

  static G4double shoot( CLHEP::HepRandomEngine* anEngine, 
                                  G4double k, G4double lambda );

  static void shootArray ( CLHEP::HepRandomEngine* anEngine, const G4int size,
                            G4double* vect, G4double k=1.0,
                            G4double lambda=1.0 );

  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  inline G4double fire();

  G4double fire( G4double k, G4double lambda );
  
  void fireArray ( const G4int size, G4double* vect);
  void fireArray ( const G4int size, G4double* vect,
                   G4double k, G4double lambda );
  inline G4double operator()();
  inline G4double operator()( G4double k, G4double lambda );

private:

  static G4double genGamma( CLHEP::HepRandomEngine *anEngine, G4double k,
                                                        G4double lambda );

  CLHEP::HepRandomEngine* localEngine;
  G4bool deleteEngine;
  G4double defaultK;
  G4double defaultLambda;

};

#include "G4MTRandGamma.icc"

#endif
