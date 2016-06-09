// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- HepStat ---
//          Purely static class containing useful statistics methods

// -----------------------------------------------------------------------

// HepStat is a substitute for using a namespace.
// One would never instantiate a HepStat object;
// usage of any of these methods looks like --
// 
// double x = HepStat::erf ( .1 );
//
// A user may wish to improve the readability of algortihm code which uses 
// one method many times by lines like using HepStat::erf
//
// and later, x = erf(u); will work.
//

// These methods are implemented in separate .cc files so that 
// user code need pull in only the code that is necessary.  Time
// (ROUGH estimates in cycles) and table footprint info is provided
// in this header.


// =======================================================================
// M. Fischler    - Created: 1/25/00
//
// M. Fischler	  - Inserted flatToGaussian 1/25/00
//			    From code of an attempt to speed up RandGauss
//			    by use of tables and splines.  The code was not
//		    	    significantly faster than Box-Mueller, so that
//		    	    algorithm is left as the RandGauss implementation.
//		  - Inserted inverseErf
// M. Fischler	  - Inserted gammln 2/4/00
// M. Fischler	  - Made constructor private; removed private destructor 4/17/00
// =======================================================================

#ifndef HepStat_h
#define HepStat_h 1

namespace CLHEP {

/**
 * @author
 * @ingroup random
 */
class HepStat {

private:
  HepStat();	
  // You CANNOT instantiate a HepStat object.

public:

  static double flatToGaussian (double r);
   // This is defined by the satement that if e() provides a uniform random
   // on (0,1) then flatToGaussian(e()) is distributed as a unit normal
   // Gaussian.  That is, flatToGaussian is the inverse of the c.d.f. of
   // a Gaussian.
  	// Footprint:  30 K  		// Time:  150 cycles

  static double inverseErf (double t);
  static double erf (double x);
        // defined in flatToGaussian.cc

  static double erfQ (double x);
  // Quicker, and with less footprint, than erf and gaussianCDF
  // but only accurate to 7 digits.
	  // Footprint:  0		// Time:  

  static double gammln (double x);
  // ln (gamma(x))

};

}  // namespace CLHEP

#endif
