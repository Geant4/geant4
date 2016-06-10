// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- RandLandau ---
//                          class header file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// Class defining methods for shooting or firing Landau distributed 
// random values.   
// 
// The Landau distribution is parameterless and describes the fluctuations
// in energy loss of a particle, making certain assumptions.  For 
// definitions and algorithms, the following papers could be read:
//
// Landau, Jour Phys VIII, No. 4, p. 201 (1944)
// Borsh-Supan, Jour Res. of NBS 65B NO. 4 p. 245 (1961)
// Kolbig & Schorr Comp Phys Comm 31 p. 97 (1984)
//
// The algorithm implemented comes form RANLAN in CERNLIB.

// =======================================================================
// M. Fischler    - Created: 5th January 2000
// M Fischler      - put and get to/from streams 12/10/04
//
// =======================================================================

#ifndef RandLandau_h
#define RandLandau_h 1

#include "CLHEP/Random/Random.h"
#include "CLHEP/Utility/memory.h"

namespace CLHEP {

/**
 * @author
 * @ingroup random
 */
class RandLandau : public HepRandom {

public:

  inline RandLandau ( HepRandomEngine& anEngine );
  inline RandLandau ( HepRandomEngine* anEngine ); 

  // These constructors should be used to instantiate a RandLandau
  // distribution object defining a local engine for it.
  // The static generator will be skipped using the non-static methods
  // defined below.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the RandLandau destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the RandLandau destructor.

  virtual ~RandLandau();
  // Destructor

  // Save and restore to/from streams
  
  std::ostream & put ( std::ostream & os ) const;
  std::istream & get ( std::istream & is );

  //
  // Methods to generate Landau-distributed random deviates.
  //
  // These deviates are accurate to the actual Landau distribution to
  // one part in 10**5 or better.

  // Static methods to shoot random values using the static generator

  static  inline double shoot();

  static  void shootArray ( const int size, double* vect );

  //  Static methods to shoot random values using a given engine
  //  by-passing the static generator.

  static  inline double shoot( HepRandomEngine* anotherEngine );

  static  void shootArray ( HepRandomEngine* anotherEngine, 
			    const int size,
                            double* vect );

  //  Instance methods using the localEngine to instead of the static 
  //  generator, and the default mean and stdDev established at construction

  inline double fire();

  void fireArray  ( const int size, double* vect);

  inline double operator()();

  std::string name() const;
  HepRandomEngine & engine();

  static std::string distributionName() {return "RandLandau";}  
  // Provides the name of this distribution class
  

protected:

  static double transform      (double r);
  static double transformSmall (double r);

private:

  std::shared_ptr<HepRandomEngine> localEngine;

};

}  // namespace CLHEP

#include "CLHEP/Random/RandLandau.icc"

#endif
