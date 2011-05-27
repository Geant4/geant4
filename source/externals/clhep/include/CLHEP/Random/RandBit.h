// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                           --- RandBit ---
//                          class header file
// -----------------------------------------------------------------------
// 

// Class defining methods for shooting Flat or Bit random numbers, double or
// integers.
// It provides methods to fill with double flat values arrays of
// specified size, as well as methods for shooting sequences of 0,1 (bits).
// Default boundaries ]0.1[ for operator()().

// This is derived from RandFlat and is a drop-in replacement.  However
// the shootBit() and fireBit() methods are stateless (which makes them
// an order of magnitude slower, but allows save/restore engine status
// to work correctly).

// =======================================================================
// M. Fischler    - Created: 15th Feb 2000
// M Fischler      - put and get to/from streams 12/10/04
// M Fischler      - static save/restore to streams streams 12/20/04
// =======================================================================

#ifndef RandBit_h
#define RandBit_h 1

#include "CLHEP/Random/RandFlat.h"

namespace CLHEP {

/**
 * @author
 * @ingroup random
 */
class RandBit : public RandFlat {

public:

  inline RandBit ( HepRandomEngine& anEngine );
  inline RandBit ( HepRandomEngine& anEngine, double width );
  inline RandBit ( HepRandomEngine& anEngine, double a, double b );
  inline RandBit ( HepRandomEngine* anEngine );
  inline RandBit ( HepRandomEngine* anEngine, double width );
  inline RandBit ( HepRandomEngine* anEngine, double a, double b );
  // These constructors should be used to instantiate a RandBit
  // distribution object defining a local engine for it.
  // The static generator will be skipped using the non-static methods
  // defined below.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the RandBit destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the RandBit destructor.

  virtual ~RandBit();
  // Destructor

  // Other than the Bit routines, constructors, and destructor, everything is
  // simply inherited from RandFlat.

  static  inline int shootBit();

  static  inline int shootBit( HepRandomEngine* );

  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  inline int fireBit();

  // Save and restore to/from streams
  
  std::ostream & put ( std::ostream & os ) const;
  std::istream & get ( std::istream & is );

  std::string name() const;

  static std::string distributionName() {return "RandBit";}  
  // Provides the name of this distribution class

  static std::ostream& saveFullState ( std::ostream & os ) 
  // Saves to stream the state of the engine and cached data.
  					{return RandFlat::saveFullState(os);}
					
  static std::istream& restoreFullState ( std::istream & is )
  // Restores from stream the state of the engine and cached data.
  					{return RandFlat::restoreFullState(is);}

  static std::ostream& saveDistState ( std::ostream & os )
  // Saves to stream the state of the cached data.
  					{return RandFlat::saveDistState(os);}

  static std::istream& restoreDistState ( std::istream & is )
  // Restores from stream the state of the cached data.
  					{return RandFlat::restoreDistState(is);}


private:

  // All the engine info, and the default A and B, are in the RandFlat
  // base class.

};

}  // namespace CLHEP

#include "CLHEP/Random/RandBit.icc"

#endif
