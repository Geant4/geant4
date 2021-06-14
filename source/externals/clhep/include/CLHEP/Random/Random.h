// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- HepRandom ---
//                          class header file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).
//
// It's a singleton instantiated by default within the HEP Random module.
// It uses an instantiated HepJamesRandom engine as default algorithm
// for pseudo-random number generation. HepRandom defines a static private
// data member theGenerator and a set of static inlined methods to manipulate
// it. By means of theGenerator the user can change the underlying engine
// algorithm, get and set the seeds and use any kind of defined random
// distribution.
// Distribution classes inherit from HepRandom and define both static and
// not-static interfaces.
// A static table of uncorrelated seeds is available in this class.
// A static method "getTheTableSeeds()" is defined to access a couple of
// seeds at a given index in the table.

// =======================================================================
// Gabriele Cosmo - Created: 5th Sep 1995
//                - Minor update: 17th May 1996
//                - Poisson now operates on doubles : 31st Oct 1996
//                - Added methods for engine status: 19th Nov 1996
//                - Fixed default values to setTheSeed() and
//                  setTheSeeds() static methods: 16th Oct 1997
//                - Modified HepRandom to act as a singleton, constructors
//                  are kept public for backward compatibility. Added table
//                  of seeds from HepRandomEngine: 19th Mar 1998
//                - Relocated Poisson and Gauss data and simplified
//                  initialisation of static generator: 5th Jan 1999
// =======================================================================

#ifndef HepRandom_h
#define HepRandom_h 1

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {

/**
 * @author <Gabriele.Cosmo@cern.ch>
 * @ingroup random
 */
class HepRandom {

public:

  HepRandom();
  HepRandom(long seed);
  // Contructors with and without a seed using the default engine
  // (MixMax).
 
  HepRandom(HepRandomEngine & algorithm);
  HepRandom(HepRandomEngine * algorithm);
  // Constructor taking an alternative engine as argument. If a pointer is
  // given the corresponding object will be deleted by the HepRandom
  // destructor.
  
  virtual ~HepRandom();
  // Destructor
  
  // implicitly allow compiler-generated copy functions 

  double flat();
  // Returns the flat value ( interval ]0...1[ ).

  void flatArray(const int size, double* vect);
  // Fills "vect" array of flat random values, given the size.

  inline double flat (HepRandomEngine* theNewEngine);
  // Returns a flat value, given a defined Random Engine.

  inline void flatArray(HepRandomEngine* theNewEngine, 
                        const int size, double* vect);
  // Fills "vect" array of flat random values, given the size
  // and a defined Random Engine.

  virtual double operator()();
  // To get a flat random number using the operator ().

  virtual std::string name() const;
  virtual HepRandomEngine & engine();
    
  
  virtual std::ostream & put ( std::ostream & os ) const;
  virtual std::istream & get ( std::istream & is );
  // Save and restore to/from streams

  // --------------------------------------------------
  // Static member functions using the static generator
  // --------------------------------------------------

  static void setTheSeed(long seed, int lxr=3);
  // (Re)Initializes the generator with a seed.

  static long getTheSeed();
  // Gets the current seed of the current generator.

  static void setTheSeeds(const long* seeds, int aux=-1);
  // (Re)Initializes the generator with a zero terminated list of seeds.

  static const long* getTheSeeds();
  // Gets the current array of seeds of the current generator.

  static void getTheTableSeeds (long* seeds, int index);
  // Gets the array of seeds in the static seedTable at "index" position.

  static HepRandom * getTheGenerator();
  // Return the current static generator.

  static void setTheEngine (HepRandomEngine* theNewEngine);
  // To set the underlying algorithm object.

  static HepRandomEngine * getTheEngine();
  // Returns a pointer to the underlying algorithm object.

  static void saveEngineStatus( const char filename[] = "Config.conf" );
  // Saves to file the current status of the current engine.

  static void restoreEngineStatus( const char filename[] = "Config.conf" );
  // Restores a saved status (if any) for the current engine.

  static std::ostream& saveFullState ( std::ostream & os );
  // Saves to stream the state of the engine and cached data.

  static std::istream& restoreFullState ( std::istream & is );
  // Restores from stream the state of the engine and cached data.

  static std::ostream& saveDistState ( std::ostream & os ) {return os;}
  // Saves to stream the state of the cached data.

  static std::istream& restoreDistState ( std::istream & is ) {return is;}
  // Restores from stream the state of the cached data.

  static std::ostream& saveStaticRandomStates ( std::ostream & os );
  // Saves to stream the engine and cached data for all distributions.

  static std::istream& restoreStaticRandomStates ( std::istream & is );
  // Restores from stream the engine and cached data for all distributions.

  static void showEngineStatus();
  // Dumps the current engine status on screen.

  static int createInstance();
  // used to initialise the default engine

  static std::string distributionName() {return "HepRandomEngine";}  
  // Provides the name of this distribution class
       
protected:     // -------- Data members ---------

  static const long seedTable[215][2];
  // Table of seeds

};

std::ostream & operator<< (std::ostream & os, const HepRandom & dist);
std::istream & operator>> (std::istream & is, HepRandom & dist);

}  // namespace CLHEP

#include "CLHEP/Random/Random.icc"

#endif
