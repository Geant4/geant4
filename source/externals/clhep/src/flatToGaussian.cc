// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- flatToGaussian ---
//                      class implementation file
// -----------------------------------------------------------------------

// Contains the methods that depend on the 30K-footpring gaussTables.cdat.
//
// flatToGaussian (double x)
// inverseErf     (double x)
// erf		  (double x)

// =======================================================================
// M Fischler	  - Created 1/25/00.
//
// =======================================================================

#include "CLHEP/Random/Stat.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include <iostream>
#include <cmath>

namespace CLHEP {

double transformSmall (double r);

//
// Table of errInts, for use with transform(r) and quickTransform(r)
//

#ifdef Traces
#define Trace1
#define Trace2
#define Trace3
#endif

// Since all these are this is static to this compilation unit only, the 
// info is establised a priori and not at each invocation.

// The main data is of course the gaussTables table; the rest is all 
// bookkeeping to know what the tables mean.

#define Table0size   200
#define Table1size   250
#define Table2size   200
#define Table3size   250
#define Table4size  1000
#define TableSize   (Table0size+Table1size+Table2size+Table3size+Table4size)

static const int Tsizes[5] =   { Table0size,
				 Table1size,
				 Table2size,
				 Table3size,
				 Table4size };

#define Table0step  (2.0E-13)
#define Table1step  (4.0E-11)  
#define Table2step  (1.0E-8) 
#define Table3step  (2.0E-6) 
#define Table4step  (5.0E-4)

static const double Tsteps[5] = { Table0step,
				 Table1step,
				 Table2step,
				 Table3step,
				 Table4step };

#define Table0offset 0
#define Table1offset (2*(Table0size))
#define Table2offset (2*(Table0size + Table1size))
#define Table3offset (2*(Table0size + Table1size + Table2size))
#define Table4offset (2*(Table0size + Table1size + Table2size + Table3size))

static const int Toffsets[5] = { Table0offset,
				 Table1offset,
				 Table2offset,
				 Table3offset,
				 Table4offset };

  // Here comes the big (30K bytes) table, kept in a file ---

static const double gaussTables [2*TableSize] = {
#include "CLHEP/Random/gaussTables.cdat"
};

double HepStat::flatToGaussian (double r) {

  double sign = +1.0;	// We always compute a negative number of 
				// sigmas.  For r > 0 we will multiply by
				// sign = -1 to return a positive number.
#ifdef Trace1
std::cout << "r = " << r << "\n";
#endif

  if ( r > .5 ) {
    r = 1-r;
    sign = -1.0;
#ifdef Trace1
std::cout << "r = " << r << "sign negative \n";
#endif
  } else if ( r == .5 ) {
    return 0.0;
  }  

  // Find a pointer to the proper table entries, along with the fraction 
  // of the way in the relevant bin dx and the bin size h.
  
  // Optimize for the case of table 4 by testing for that first.  
  // By removing that case from the for loop below, we save not only
  // several table lookups, but also several index calculations that
  // now become (compile-time) constants.
  //
  // Past the case of table 4, we need not be as concerned about speed since
  // this will happen only .1% of the time.

  const double* tptr = 0;
  double  dx = 0;
  double  h = 0;

  // The following big if block will locate tptr, h, and dx from whichever
  // table is applicable:

  int index;

  if ( r >= Table4step ) {

    index = int((Table4size<<1) * r);	// 1 to Table4size-1 
    if (index <= 0) index = 1; 			// in case of rounding problem
    if (index >= Table4size) index = Table4size-1;
    dx = (Table4size<<1) * r - index; 		// fraction of way to next bin
    h = Table4step;
#ifdef Trace2 
std::cout << "index = " << index << " dx = " << dx << " h = " << h << "\n";
#endif
    index = (index<<1) + (Table4offset-2);	
	// at r = table4step+eps, index refers to the start of table 4 
	// and at r = .5 - eps, index refers to the next-to-last entry:
	// remember, the table has two numbers per actual entry.
#ifdef Trace2 
std::cout << "offset index = " << index << "\n";
#endif

    tptr = &gaussTables [index];
    
  } else if (r < Tsteps[0])  {

    // If r is so small none of the tables apply, use the asymptotic formula.
    return (sign * transformSmall (r));

  } else {
    
    for ( int tableN = 3; tableN >= 0; tableN-- ) {
      if ( r < Tsteps[tableN] ) {continue;} 	// can't happen when tableN==0
#ifdef Trace2 
std::cout << "Using table " << tableN << "\n";
#endif
      double step = Tsteps[tableN];
      index = int(r/step);			// 1 to TableNsize-1 
        // The following two tests should probably never be true, but
        // roundoff might cause index to be outside its proper range.
        // In such a case, the interpolation still makes sense, but we
        // need to take care that tptr points to valid data in the right table.
      if (index == 0) index = 1; 			
      if (index >= Tsizes[tableN]) index = Tsizes[tableN] - 1;
      dx =  r/step - index; 			// fraction of way to next bin
      h  =  step;
#ifdef Trace2 
std::cout << "index = " << index << " dx = " << dx << " h = " << h << "\n";
#endif
      index = (index<<1) + Toffsets[tableN] - 2;
      tptr = &gaussTables [index];
      break;
    } // end of the for loop which finds tptr, dx and h in one of the tables

    // The code can only get to here by a break statement, having set dx etc.
    // It not get to here without going through one of the breaks,
    // because before the for loop we test for the case of r < Tsteps[0].

  } // End of the big if block.

  // At the end of this if block, we have tptr, dx and h -- and if r is less 
  // than the smallest step, we have already returned the proper answer.  

  double  y0 = *tptr++;
  double  d0 = *tptr++;
  double  y1 = *tptr++;
  double  d1 = *tptr;

#ifdef Trace3
std::cout << "y0: " << y0 << " y1: " << y1 << " d0: " << d0 << " d1: " << d1 << "\n";
#endif

  double  x2 = dx * dx;
  double  oneMinusX = 1 - dx;
  double  oneMinusX2 = oneMinusX * oneMinusX;

  double  f0 = (2. * dx + 1.) * oneMinusX2;
  double  f1 = (3. - 2. * dx) * x2;
  double  g0 =  h * dx * oneMinusX2;
  double  g1 =  - h * oneMinusX * x2;

#ifdef Trace3
std::cout << "f0: " << f0 << " f1: " << f1 << " g0: " << g0 << " g1: " << g1 << "\n";
#endif

  double answer = f0 * y0 + f1 * y1 + g0 * d0 + g1 * d1;

#ifdef Trace1
std::cout << "variate is: " << sign*answer << "\n";
#endif

  return sign * answer;

} // flatToGaussian

double transformSmall (double r) {

  // Solve for -v in the asymtotic formula 
  //
  // errInt (-v) =  exp(-v*v/2)         1     1*3    1*3*5
  //		   ------------ * (1 - ---- + ---- - ----- + ... )
  //		   v*sqrt(2*pi)        v**2   v**4   v**6

  // The value of r (=errInt(-v)) supplied is going to less than 2.0E-13,
  // which is such that v < -7.25.  Since the value of r is meaningful only
  // to an absolute error of 1E-16 (double precision accuracy for a number 
  // which on the high side could be of the form 1-epsilon), computing
  // v to more than 3-4 digits of accuracy is suspect; however, to ensure 
  // smoothness with the table generator (which uses quite a few terms) we
  // also use terms up to 1*3*5* ... *13/v**14, and insist on accuracy of
  // solution at the level of 1.0e-7.

  // This routine is called less than one time in a trillion firings, so
  // speed is of no concern.  As a matter of technique, we terminate the
  // iterations in case they would be infinite, but this should not happen.

  double eps = 1.0e-7;
  double guess = 7.5;
  double v;
  
  for ( int i = 1; i < 50; i++ ) {
    double vn2 = 1.0/(guess*guess);
    double s1 = -13*11*9*7*5*3 * vn2*vn2*vn2*vn2*vn2*vn2*vn2;
            s1 +=    11*9*7*5*3 * vn2*vn2*vn2*vn2*vn2*vn2;
            s1 +=      -9*7*5*3 * vn2*vn2*vn2*vn2*vn2;
            s1 +=         7*5*3 * vn2*vn2*vn2*vn2;
            s1 +=          -5*3 * vn2*vn2*vn2;
            s1 +=            3 * vn2*vn2    - vn2  +    1.0;
    v = std::sqrt ( 2.0 * std::log ( s1 / (r*guess*std::sqrt(CLHEP::twopi)) ) );
    if ( std::abs(v-guess) < eps ) break;
    guess = v;
  }
 
  return -v;

} // transformSmall()

double HepStat::inverseErf (double t) {

  // This uses erf(x) = 2*gaussCDF(sqrt(2)*x) - 1

  return std::sqrt(0.5) * flatToGaussian(.5*(t+1));

}

double HepStat::erf (double x) {

// Math:
//
// For any given x we can "quickly" find t0 = erfQ (x) = erf(x) + epsilon.
//
// Then we can find x1 = inverseErf (t0) = inverseErf (erf(x)+epsilon)
//
// Expanding in the small epsion, 
// 
//  x1 = x + epsilon * [deriv(inverseErf(u),u) at u = t0] + O(epsilon**2)
//
// so epsilon is (x1-x) / [deriv(inverseErf(u),u) at u = t0] + O(epsilon**2)
//
// But the derivative of an inverse function is one over the derivative of the
// function, so 
// epsilon  = (x1-x) * [deriv(erf(v),v) at v = x] + O(epsilon**2)
//
// And the definition of the erf integral makes that derivative trivial.
// Ultimately,
//
// erf(x) = erfQ(x) - (inverseErf(erfQ(x))-x) * exp(-x**2) * 2/sqrt(CLHEP::pi)
//

  double t0 = erfQ(x);
  double deriv = std::exp(-x*x) * (2.0 / std::sqrt(CLHEP::pi));

  return t0 - (inverseErf (t0) - x) * deriv;

}


}  // namespace CLHEP

