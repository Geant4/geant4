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
// Class G4MTRandFlat
// 
// Modified version with MT extensions
// of corresponding CLHEP class RandFlat

// --------------------------------------------------------------------
#ifndef G4MTRandFlat_hh
#define G4MTRandFlat_hh 1

#include "G4MTHepRandom.hh"

class G4MTRandFlat : public G4MTHepRandom
{

public:

  inline G4MTRandFlat ( CLHEP::HepRandomEngine& anEngine );
  inline G4MTRandFlat ( CLHEP::HepRandomEngine& anEngine, G4double width );
  inline G4MTRandFlat ( CLHEP::HepRandomEngine& anEngine, G4double a, G4double b );
  inline G4MTRandFlat ( CLHEP::HepRandomEngine* anEngine );
  inline G4MTRandFlat ( CLHEP::HepRandomEngine* anEngine, G4double width );
  inline G4MTRandFlat ( CLHEP::HepRandomEngine* anEngine, G4double a, G4double b );
  // These constructors should be used to instantiate a RandFlat
  // distribution object defining a local engine for it.
  // The static generator will be skipped using the non-static methods
  // defined below.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the RandFlat destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the RandFlat destructor.

  virtual ~G4MTRandFlat();
  // Destructor

  // Static methods to shoot random values using the static generator

  static  G4double shoot();

  static  inline G4double shoot( G4double width );

  static  inline G4double shoot( G4double a, G4double b );

  static  inline G4long shootInt( G4long n );

  static  inline G4long shootInt( G4long m, G4long n );

  static  inline G4int shootBit();

  static  void shootArray ( const G4int size, G4double* vect );

  static  void shootArray ( const G4int size, G4double* vect,
                            G4double lx, G4double dx );

  //  Static methods to shoot random values using a given engine
  //  by-passing the static generator.

  static  inline G4double shoot ( CLHEP::HepRandomEngine* anEngine );

  static  inline G4double shoot( CLHEP::HepRandomEngine* anEngine, G4double width );

  static  inline G4double shoot( CLHEP::HepRandomEngine* anEngine,
                                  G4double a, G4double b );
  static  inline G4long shootInt( CLHEP::HepRandomEngine* anEngine, G4long n );
  
  static  inline G4long shootInt( CLHEP::HepRandomEngine* anEngine, G4long m, G4long n );
  
  static  inline G4int shootBit( CLHEP::HepRandomEngine* );

  static  inline void shootArray ( CLHEP::HepRandomEngine* anEngine,
                                   const G4int size, G4double* vect );

  static  void shootArray ( CLHEP::HepRandomEngine* anEngine, 
                            const G4int size, G4double* vect,
                            G4double lx, G4double dx );

  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  inline G4double fire();

  inline G4double fire( G4double width );

  inline G4double fire( G4double a, G4double b );

  inline G4long fireInt( G4long n );

  inline G4long fireInt( G4long m, G4long n );

  inline G4int fireBit();

  void fireArray (const G4int size, G4double* vect);

  void fireArray (const G4int size, G4double* vect,
                  G4double lx, G4double dx);

  G4double operator()();
  G4double operator()( G4double width );
  G4double operator()( G4double a, G4double b );

private:

  // ShootBits generates an integer random number,
  // which is used by fireBit().
  // The number is stored in randomInt and firstUnusedBit

  inline void fireBits();
  static inline void shootBits();
  static inline void shootBits(CLHEP::HepRandomEngine*);

  // In MSB, the most significant bit of the integer random number
  // generated by ShootBits() is set.
  // Note:
  //   the number of significant bits must be chosen so that
  //   - an unsigned long can hold it
  //   - and it should be less than the number of bits returned 
  //     by Shoot() which are not affected by precision problems
  //     on _each_ architecture.
  //   (Aim: the random generators should be machine-independent).

  static const unsigned long MSB; 
  static const G4int MSBBits;
  // These two are set up in RandFlat.cc and need not be saved/restored

  unsigned long randomInt;
  unsigned long firstUnusedBit;
  static G4ThreadLocal unsigned long staticRandomInt;
  static G4ThreadLocal unsigned long staticFirstUnusedBit;
  
  CLHEP::HepRandomEngine* localEngine;
  G4bool deleteEngine;
  G4double defaultA;
  G4double defaultB;

};

#include "G4MTRandFlat.icc"

#endif
