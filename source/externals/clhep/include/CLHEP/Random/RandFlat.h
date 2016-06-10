// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                           --- RandFlat ---
//                          class header file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// Class defining methods for shooting flat random numbers, double or
// integers.
// It provides methods to fill with double flat values arrays of
// specified size, as well as methods for shooting sequences of 0,1 (bits).
// Default boundaries ]0.1[ for operator()().

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
// Peter Urban    - ShootBit() and related stuff added: 5th Sep 1996
// Gabriele Cosmo - Added operator() and additional methods to fill
//                  arrays specifying boundaries: 24th Jul 1997 
// J.Marraffino   - Added default arguments as attributes and
//                  operator() with arguments: 16th Feb 1998
// M. Fischler    - Moved copy constructor to protected so that
//		    derived RandBit can get at it.
// M Fischler      - put and get to/from streams 12/10/04
// =======================================================================

#ifndef RandFlat_h
#define RandFlat_h 1

#include "CLHEP/Random/Random.h"
#include "CLHEP/Utility/memory.h"
#include "CLHEP/Utility/thread_local.h"

namespace CLHEP {

/**
 * @author <Gabriele.Cosmo@cern.ch>
 * @ingroup random
 */
class RandFlat : public HepRandom {

public:

  inline RandFlat ( HepRandomEngine& anEngine );
  inline RandFlat ( HepRandomEngine& anEngine, double width );
  inline RandFlat ( HepRandomEngine& anEngine, double a, double b );
  inline RandFlat ( HepRandomEngine* anEngine );
  inline RandFlat ( HepRandomEngine* anEngine, double width );
  inline RandFlat ( HepRandomEngine* anEngine, double a, double b );
  // These constructors should be used to instantiate a RandFlat
  // distribution object defining a local engine for it.
  // The static generator will be skipped using the non-static methods
  // defined below.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the RandFlat destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the RandFlat destructor.

  virtual ~RandFlat();
  // Destructor

  // Static methods to shoot random values using the static generator

  static  double shoot();

  static  inline double shoot( double width );

  static  inline double shoot( double a, double b );

  static  inline long shootInt( long n );

  static  inline long shootInt( long a1, long n );

  static  inline int shootBit();

  static  void shootArray ( const int size, double* vect );

  static  void shootArray ( const int size, double* vect,
                            double lx, double dx );

  //  Static methods to shoot random values using a given engine
  //  by-passing the static generator.

  static  inline double shoot ( HepRandomEngine* anEngine );

  static  inline double shoot( HepRandomEngine* anEngine, double width );

  static  inline double shoot( HepRandomEngine* anEngine,
                                  double a, double b );
  static  inline long shootInt( HepRandomEngine* anEngine, long n );
  
  static  inline long shootInt( HepRandomEngine* anEngine, long a1, long n );
  
  static  inline int shootBit( HepRandomEngine* );

  static  inline void shootArray ( HepRandomEngine* anEngine,
                                   const int size, double* vect );

  static  void shootArray ( HepRandomEngine* anEngine, 
                            const int size, double* vect,
                            double lx, double dx );

  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  inline double fire();

  inline double fire( double width );

  inline double fire( double a, double b );

  inline long fireInt( long n );

  inline long fireInt( long a1, long n );

  inline int fireBit();

  void fireArray (const int size, double* vect);

  void fireArray (const int size, double* vect,
                  double lx, double dx);

  double operator()();
  double operator()( double width );
  double operator()( double a, double b );

  // Save and restore to/from streams
  
  std::ostream & put ( std::ostream & os ) const;
  std::istream & get ( std::istream & is );

  std::string name() const;
  HepRandomEngine & engine();

  static std::string distributionName() {return "RandFlat";}  
  // Provides the name of this distribution class 
  
  // Methods overriding the base class static saveEngineStatus ones,
  // by adding extra data so that save in one program, then further shootBit()s
  // will produce the identical sequence to restore in another program, then
  // generating shootBit() randoms there

  static void saveEngineStatus( const char filename[] = "Config.conf" );
  // Saves to file the current status of the current engine.

  static void restoreEngineStatus( const char filename[] = "Config.conf" );
  // Restores a saved status (if any) for the current engine.

  static std::ostream& saveFullState ( std::ostream & os );
  // Saves to stream the state of the engine and cached data.

  static std::istream& restoreFullState ( std::istream & is );
  // Restores from stream the state of the engine and cached data.

  static std::ostream& saveDistState ( std::ostream & os );
  // Saves to stream the state of the cached data.

  static std::istream& restoreDistState ( std::istream & is );
  // Restores from stream the state of the cached data.


protected:

#if 0
  // Protected copy constructor. Defining it here disallows use by users.
  RandFlat(const RandFlat& d);
#endif  // 0

private:

  // ShootBits generates an integer random number,
  // which is used by fireBit().
  // The number is stored in randomInt and firstUnusedBit

  inline void fireBits();
  static inline void shootBits();
  static inline void shootBits(HepRandomEngine*);

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
  static const int MSBBits;
  // These two are set up in RandFlat.cc and need not be saved/restored

  unsigned long randomInt;
  unsigned long firstUnusedBit;
  static CLHEP_THREAD_LOCAL unsigned long staticRandomInt;
  static CLHEP_THREAD_LOCAL unsigned long staticFirstUnusedBit;
  
  std::shared_ptr<HepRandomEngine> localEngine;
  double defaultWidth;
  double defaultA;
  double defaultB;

};

}  // namespace CLHEP

#include "CLHEP/Random/RandFlat.icc"

#endif
