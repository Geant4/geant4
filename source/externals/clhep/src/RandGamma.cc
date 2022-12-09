// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- RandGamma ---
//                      class implementation file
// -----------------------------------------------------------------------

// =======================================================================
// John Marraffino - Created: 12th May 1998
// M Fischler      - put and get to/from streams 12/13/04
// M Fischler	      - put/get to/from streams uses pairs of ulongs when
//			+ storing doubles avoid problems with precision 
//			4/14/05
// =======================================================================

#include "CLHEP/Random/RandGamma.h"
#include "CLHEP/Random/DoubConv.h"
#include "CLHEP/Utility/thread_local.h"
#include <cmath>	// for std::log()
#include <iostream>
#include <string>
#include <vector>

namespace CLHEP {

std::string RandGamma::name() const {return "RandGamma";}
HepRandomEngine & RandGamma::engine() {return *localEngine;}

RandGamma::~RandGamma() {
}

double RandGamma::shoot( HepRandomEngine *anEngine,  double k,
                                                        double lambda ) {
  return genGamma( anEngine, k, lambda );
}

double RandGamma::shoot( double k, double lambda ) {
  HepRandomEngine *anEngine = HepRandom::getTheEngine();
  return genGamma( anEngine, k, lambda );
}

double RandGamma::fire( double k, double lambda ) {
  return genGamma( localEngine.get(), k, lambda );
}

void RandGamma::shootArray( const int size, double* vect,
                            double k, double lambda )
{
  for( double* v = vect; v != vect + size; ++v )
    *v = shoot(k,lambda);
}

void RandGamma::shootArray( HepRandomEngine* anEngine,
                            const int size, double* vect,
                            double k, double lambda )
{
  for( double* v = vect; v != vect + size; ++v )
    *v = shoot(anEngine,k,lambda);
}

void RandGamma::fireArray( const int size, double* vect)
{
  for( double* v = vect; v != vect + size; ++v )
    *v = fire(defaultK,defaultLambda);
}

void RandGamma::fireArray( const int size, double* vect,
                           double k, double lambda )
{
  for( double* v = vect; v != vect + size; ++v )
    *v = fire(k,lambda);
}

double RandGamma::genGamma( HepRandomEngine *anEngine,
                               double a, double lambda ) {
/*************************************************************************
 *         Gamma Distribution - Rejection algorithm gs combined with     *
 *                              Acceptance complement method gd          *
 *************************************************************************/

  double aa = -1.0, aaa = -1.0, b{0.}, c{0.}, d{0.}, e{0.}, r{0.}, s{0.}, si{0.}, ss{0.}, q0{0.};
  constexpr double q1 = 0.0416666664, q2 =  0.0208333723, q3 = 0.0079849875,
       q4 = 0.0015746717, q5 = -0.0003349403, q6 = 0.0003340332,
       q7 = 0.0006053049, q8 = -0.0004701849, q9 = 0.0001710320,
       a1 = 0.333333333,  a2 = -0.249999949,  a3 = 0.199999867,
       a4 =-0.166677482,  a5 =  0.142873973,  a6 =-0.124385581,
       a7 = 0.110368310,  a8 = -0.112750886,  a9 = 0.104089866,
       e1 = 1.000000000,  e2 =  0.499999994,  e3 = 0.166666848,
       e4 = 0.041664508,  e5 =  0.008345522,  e6 = 0.001353826,
       e7 = 0.000247453;

 double gds{0.},p{0.},q{0.},t{0.},sign_u{0.},u{0.},v{0.},w{0.},x{0.};
 double v1{0.},v2{0.},v12{0.};

// Check for invalid input values

 if( a <= 0.0 ) return (-1.0);
 if( lambda <= 0.0 ) return (-1.0);

 if (a < 1.0)
   {          // CASE A: Acceptance rejection algorithm gs
    b = 1.0 + 0.36788794412 * a;       // Step 1
    for(;;)
      {
       p = b * anEngine->flat();
       if (p <= 1.0)
	  {                            // Step 2. Case gds <= 1
	   gds = std::exp(std::log(p) / a);
	   if (std::log(anEngine->flat()) <= -gds) return(gds/lambda);
	  }
       else
	  {                            // Step 3. Case gds > 1
	   gds = - std::log ((b - p) / a);
	   if (std::log(anEngine->flat()) <= ((a - 1.0) * std::log(gds))) return(gds/lambda);
	  }
      }
   }
 else
   {          // CASE B: Acceptance complement algorithm gd
    if (a != aa)
       {                               // Step 1. Preparations
	aa = a;
	ss = a - 0.5;
	s = std::sqrt(ss);
	d = 5.656854249 - 12.0 * s;
       }
                                              // Step 2. Normal deviate
    do {
      v1 = 2.0 * anEngine->flat() - 1.0;
      v2 = 2.0 * anEngine->flat() - 1.0;
      v12 = v1*v1 + v2*v2;
    } while ( v12 > 1.0 );
    t = v1*std::sqrt(-2.0*std::log(v12)/v12);
    x = s + 0.5 * t;
    gds = x * x;
    if (t >= 0.0) return(gds/lambda);         // Immediate acceptance

    u = anEngine->flat();            // Step 3. Uniform random number
    if (d * u <= t * t * t) return(gds/lambda); // Squeeze acceptance

    if (a != aaa)
       {                               // Step 4. Set-up for hat case
	aaa = a;
	r = 1.0 / a;
	q0 = ((((((((q9 * r + q8) * r + q7) * r + q6) * r + q5) * r + q4) *
			  r + q3) * r + q2) * r + q1) * r;
	if (a > 3.686)
	   {
	    if (a > 13.022)
	       {
		b = 1.77;
		si = 0.75;
		c = 0.1515 / s;
	       }
	    else
	       {
		b = 1.654 + 0.0076 * ss;
		si = 1.68 / s + 0.275;
		c = 0.062 / s + 0.024;
	       }
	   }
	else
	   {
	    b = 0.463 + s - 0.178 * ss;
	    si = 1.235;
	    c = 0.195 / s - 0.079 + 0.016 * s;
	   }
       }
    if (x > 0.0)                       // Step 5. Calculation of q
       {
	v = t / (s + s);               // Step 6.
	if (std::fabs(v) > 0.25)
	   {
	    q = q0 - s * t + 0.25 * t * t + (ss + ss) * std::log(1.0 + v);
	   }
	else
	   {
	    q = q0 + 0.5 * t * t * ((((((((a9 * v + a8) * v + a7) * v + a6) *
			    v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
	   }                // Step 7. Quotient acceptance
	if (std::log(1.0 - u) <= q) return(gds/lambda);
       }

    for(;;)
       {                    // Step 8. Double exponential deviate t
	do
	{
	 e = -std::log(anEngine->flat());
	 u = anEngine->flat();
	 u = u + u - 1.0;
	 sign_u = (u > 0)? 1.0 : -1.0;
	 t = b + (e * si) * sign_u;
	}
	while (t <= -0.71874483771719);   // Step 9. Rejection of t
	v = t / (s + s);                  // Step 10. New q(t)
	if (std::fabs(v) > 0.25)
	   {
	    q = q0 - s * t + 0.25 * t * t + (ss + ss) * std::log(1.0 + v);
	   }
	else
	   {
	    q = q0 + 0.5 * t * t * ((((((((a9 * v + a8) * v + a7) * v + a6) *
			    v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
	   }
	if (q <= 0.0) continue;           // Step 11.
	if (q > 0.5)
	   {
	    w = std::exp(q) - 1.0;
	   }
	else
	   {
	    w = ((((((e7 * q + e6) * q + e5) * q + e4) * q + e3) * q + e2) *
				     q + e1) * q;
	   }                    // Step 12. Hat acceptance
	if ( c * u * sign_u <= w * std::exp(e - 0.5 * t * t))
	   {
	    x = s + 0.5 * t;
	    return(x*x/lambda);
	   }
       }
   }
}

std::ostream & RandGamma::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  t = DoubConv::dto2longs(defaultK);
  os << defaultK << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(defaultLambda);
  os << defaultLambda << " " << t[0] << " " << t[1] << "\n";
  os.precision(pr);
  return os;
}

std::istream & RandGamma::get ( std::istream & is ) {
  std::string inName;
  is >> inName;
  if (inName != name()) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Mismatch when expecting to read state of a "
    	      << name() << " distribution\n"
	      << "Name found was " << inName
	      << "\nistream is left in the badbit state\n";
    return is;
  }
  if (possibleKeywordInput(is, "Uvec", defaultK)) {
    std::vector<unsigned long> t(2);
    is >> defaultK >> t[0] >> t[1]; defaultK = DoubConv::longs2double(t); 
    is >> defaultLambda>>t[0]>>t[1]; defaultLambda = DoubConv::longs2double(t); 
    return is;
  }
  // is >> defaultK encompassed by possibleKeywordInput
  is >> defaultLambda;
  return is;
}

}  // namespace CLHEP

