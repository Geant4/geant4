// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                         --- RandStudentT ---
//                          class header file
// -----------------------------------------------------------------------

// Class defining methods for shooting Student's t- distributed random 
// values, given a number of degrees of freedom a (default=1.0).
// Default values are used for operator()().

// Valid input values are a > 0.  When invalid values are presented, the
// code silently returns DBL_MAX from <float.h> which is the same as
// MAXDOUBLE in <values.h> on systems where the latter exists.

// =======================================================================
// John Marraffino - Created: Based on the C-Rand package
//                   by Ernst Stadlober and Franz Niederl of the Technical
//                   University of Graz, Austria : 12th May 1998
//                 - Removed <values.h> because that won't work
//                   on NT : 26th Jun 1998
// Gabriele Cosmo  - Fixed minor bug on inline definition for shoot()
//                   methods. Created .icc file : 20th Aug 1998
//                 - Removed useless methods and data: 5th Jan 1999
// M Fischler      - put and get to/from streams 12/10/04
// =======================================================================

#ifndef RandStudentT_h
#define RandStudentT_h 1

#include "CLHEP/Random/Random.h"
#include "CLHEP/Utility/memory.h"

namespace CLHEP {

/**
 * @author
 * @ingroup random
 */
class RandStudentT : public HepRandom {

public:

  inline RandStudentT ( HepRandomEngine& anEngine, double a=1.0 );
  inline RandStudentT ( HepRandomEngine* anEngine, double a=1.0 );
  // These constructors should be used to instantiate a RandStudentT
  // distribution object defining a local engine for it.
  // The static generator will be skipped using the non-static methods
  // defined below.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the RandStudentT destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the RandStudentT destructor.

  virtual ~RandStudentT();
  // Destructor

  // Save and restore to/from streams
  
  std::ostream & put ( std::ostream & os ) const;
  std::istream & get ( std::istream & is );

  // Static methods to shoot random values using the static generator

  static  inline double shoot();

  static  double shoot( double a );

  static  void shootArray ( const int size, double* vect,
                            double a=1.0 );

  //  Static methods to shoot random values using a given engine
  //  by-passing the static generator.

  static  inline double shoot( HepRandomEngine* anEngine );

  static  double shoot( HepRandomEngine* anEngine, 
                           double a );

  static  void shootArray ( HepRandomEngine* anEngine, const int size,
                            double* vect, double a=1.0 );

  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  inline double fire();

  double fire( double a );
  
  void fireArray ( const int size, double* vect );
  void fireArray ( const int size, double* vect, double a );
  double operator()();
  double operator()( double a );

  std::string name() const;
  HepRandomEngine & engine();

  static std::string distributionName() {return "RandStudentT";}  
  // Provides the name of this distribution class


private:

  std::shared_ptr<HepRandomEngine> localEngine;
  double defaultA;

};

}  // namespace CLHEP

#include "CLHEP/Random/RandStudentT.icc"

#endif
