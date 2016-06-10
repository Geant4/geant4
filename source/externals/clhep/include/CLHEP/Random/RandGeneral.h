// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- RandGeneral ---
//                          class header file
// -----------------------------------------------------------------------

// Class defining methods for shooting generally distributed random values,
// given a user-defined probability distribution function.

// =======================================================================
// S.Magni & G.Pieri  - Created: 29 April 1998 
// G.Cosmo            - Added constructor using default engine from the
//                     	static generator: 20 Aug 1998
// S.Magni & G.Pieri  - Added linear interpolation: 24 March 1999
// M. Fischler	      - Added private methods that simplify the implementaion
// 			prepareTables(), useFlatDistribution(), mapRandom()
//		      - Added private variable oneOverNbins.
//	     	      - Made the warning about shoot() not being static a tad
//			more prominent.   		14 May 1999	
// M Fischler         - put and get to/from streams 12/15/04
// =======================================================================

#ifndef RandGeneral_h
#define RandGeneral_h 1

#include "CLHEP/Random/Random.h"
#include "CLHEP/Utility/memory.h"
#include <vector>

namespace CLHEP {

/**
 * @author
 * @ingroup random
 */
class RandGeneral : public HepRandom {

public:

  RandGeneral ( const double* aProbFunc, 
		int theProbSize, 
		int IntType=0 );
  RandGeneral ( HepRandomEngine& anEngine,
                const double* aProbFunc, 
		int theProbSize, 
		int IntType=0 );
  RandGeneral ( HepRandomEngine* anEngine, 
                const double* aProbFunc, 
		int theProbSize, 
		int IntType=0 );
  // These constructors should be used to instantiate a RandGeneral
  // distribution object defining a local engine for it.
  // The static generator will be skipped by using the non-static methods
  // defined below. In case no engine is specified in the constructor, the
  // default engine used by the static generator is applied.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the RandGeneral destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the RandGeneral destructor.
  // The probability distribution function (Pdf) must be provided by the user
  // as an array of positive real number. The array size must also be
  // provided. The Pdf doesn't need to be normalized to 1. 
  // if IntType = 0 ( default value ) a uniform random number is
  // generated using the engine. The uniform number is then transformed
  // to the user's distribution using the cumulative probability
  // distribution constructed from his histogram. The cumulative
  // distribution is inverted using a binary search for the nearest
  // bin boundary and a linear interpolation within the
  // bin. RandGeneral therefore generates a constant density within
  // each bin.
  // if IntType = 1 no interpolation is performed and the result is a
  // discrete distribution.

  virtual ~RandGeneral();
  // Destructor

  // Methods to shoot random values using the static generator
  // N.B.: The methods are NOT static since they use nonstatic members
  // theIntegralPdf & nBins

	/////////////////////
	//		   //
	// BIG RED WARNING //
	//		   //
	/////////////////////
	//
	// The above N.B. is telling users that the shoot() methods in this
	// class are NOT STATIC.  You cannot do 
	//	double x = RandGeneral::shoot();
	// It would not make sense to provide a static shoot -- what would 
	// the default probability function look like?

  inline double shoot();

  inline void shootArray ( const int size, double* vect);

  //  Methods to shoot random values using a given engine
  //  by-passing the static generator.

  double shoot( HepRandomEngine* anEngine );

  void shootArray ( HepRandomEngine* anEngine, const int size,
                    double* vect );
			    
  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  double fire();

  void fireArray ( const int size, double* vect);

  double operator()();

  // Save and restore to/from streams
  
  std::ostream & put ( std::ostream & os ) const;
  std::istream & get ( std::istream & is );

  std::string name() const;
  HepRandomEngine & engine();

  static std::string distributionName() {return "RandGeneral";}  
  // Provides the name of this distribution class
  

private:

  std::shared_ptr<HepRandomEngine> localEngine;
  std::vector<double> theIntegralPdf;
  int nBins;
  double oneOverNbins;
  int InterpolationType;

  // Private methods to factor out replicated implementation sections
  void prepareTable(const double* aProbFunc);
  void useFlatDistribution();
  double mapRandom(double rand) const;

};

}  // namespace CLHEP

#include "CLHEP/Random/RandGeneral.icc"

#endif
