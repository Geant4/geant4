// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- RandChiSquare ---
//                      class implementation file
// -----------------------------------------------------------------------

// =======================================================================
// John Marraffino - Created: 12th May 1998
// M Fischler     - put and get to/from streams 12/10/04
// M Fischler	      - put/get to/from streams uses pairs of ulongs when
//			+ storing doubles avoid problems with precision 
//			4/14/05
// =======================================================================

#include "CLHEP/Random/RandChiSquare.h"
#include "CLHEP/Random/DoubConv.h"
#include "CLHEP/Utility/thread_local.h"
#include <cmath>	// for std::log()

namespace CLHEP {

std::string RandChiSquare::name() const {return "RandChiSquare";}
HepRandomEngine & RandChiSquare::engine() {return *localEngine;}

RandChiSquare::~RandChiSquare() {
}

double RandChiSquare::shoot( HepRandomEngine *anEngine,  double a ) {
  return genChiSquare( anEngine, a );
}

double RandChiSquare::shoot( double a ) {
  HepRandomEngine *anEngine = HepRandom::getTheEngine();
  return genChiSquare( anEngine, a );
}

double RandChiSquare::fire( double a ) {
  return genChiSquare( localEngine.get(), a );
}

void RandChiSquare::shootArray( const int size, double* vect,
                            double a ) {
  for( double* v = vect; v != vect+size; ++v )
    *v = shoot(a);
}

void RandChiSquare::shootArray( HepRandomEngine* anEngine,
                            const int size, double* vect,
                            double a )
{
  for( double* v = vect; v != vect+size; ++v )
    *v = shoot(anEngine,a);
}

void RandChiSquare::fireArray( const int size, double* vect) {
  for( double* v = vect; v != vect+size; ++v )
    *v = fire(defaultA);
}

void RandChiSquare::fireArray( const int size, double* vect,
                           double a ) {
  for( double* v = vect; v != vect+size; ++v )
    *v = fire(a);
}

double RandChiSquare::genChiSquare( HepRandomEngine *anEngine,
                                       double a ) {
/******************************************************************
 *                                                                *
 *        Chi Distribution - Ratio of Uniforms  with shift        *
 *                                                                *
 ******************************************************************
 *                                                                *
 * FUNCTION :   - chru samples a random number from the Chi       *
 *                distribution with parameter  a > 1.             *
 * REFERENCE :  - J.F. Monahan (1987): An algorithm for           *
 *                generating chi random variables, ACM Trans.     *
 *                Math. Software 13, 168-172.                     *
 * SUBPROGRAM : - anEngine  ... pointer to a (0,1)-Uniform        *
 *                engine                                          *
 *                                                                *
 * Implemented by R. Kremer, 1990                                 *
 ******************************************************************/

 static CLHEP_THREAD_LOCAL double a_in = -1.0,b,vm,vp,vd;
 double u,v,z,zz,r;

// Check for invalid input value

 if( a < 1 )  return (-1.0);

 if (a == 1)
  {
   for(;;)
    {
     u = anEngine->flat();
     v = anEngine->flat() * 0.857763884960707;
     z = v / u;
     if (z < 0) continue;
     zz = z * z;
     r = 2.5 - zz;
     if (z < 0.0) r = r + zz * z / (3.0 * z);
     if (u < r * 0.3894003915) return(z*z);
     if (zz > (1.036961043 / u + 1.4)) continue;
     if (2 * std::log(u) < (- zz * 0.5 )) return(z*z);
     }
   }
 else
  {
   if (a != a_in)
    {
     b = std::sqrt(a - 1.0);
     vm = - 0.6065306597 * (1.0 - 0.25 / (b * b + 1.0));
     vm = (-b > vm)? -b : vm;
     vp = 0.6065306597 * (0.7071067812 + b) / (0.5 + b);
     vd = vp - vm;
     a_in = a;
     }
   for(;;)
    {
     u = anEngine->flat();
     v = anEngine->flat() * vd + vm;
     z = v / u;
     if (z < -b) continue;
     zz = z * z;
     r = 2.5 - zz;
     if (z < 0.0) r = r + zz * z / (3.0 * (z + b));
     if (u < r * 0.3894003915) return((z + b)*(z + b));
     if (zz > (1.036961043 / u + 1.4)) continue;
     if (2 * std::log(u) < (std::log(1.0 + z / b) * b * b - zz * 0.5 - z * b)) return((z + b)*(z + b));
     }
   }
}

std::ostream & RandChiSquare::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  t = DoubConv::dto2longs(defaultA);
  os << defaultA << " " << t[0] << " " << t[1] << "\n";
  os.precision(pr);
  return os;
}

std::istream & RandChiSquare::get ( std::istream & is ) {
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
  if (possibleKeywordInput(is, "Uvec", defaultA)) {
    std::vector<unsigned long> t(2);
    is >> defaultA >> t[0] >> t[1]; defaultA = DoubConv::longs2double(t); 
    return is;
  }
  // is >> defaultA encompassed by possibleKeywordInput
  return is;
}



}  // namespace CLHEP

