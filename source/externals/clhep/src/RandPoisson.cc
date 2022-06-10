// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                         --- RandPoisson ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
//                - Added not static Shoot() method: 17th May 1996
//                - Algorithm now operates on doubles: 31st Oct 1996
//                - Added methods to shoot arrays: 28th July 1997
//                - Added check in case xm=-1: 4th February 1998
// J.Marraffino   - Added default mean as attribute and
//                  operator() with mean: 16th Feb 1998
// Gabriele Cosmo - Relocated static data from HepRandom: 5th Jan 1999
// M Fischler     - put and get to/from streams 12/15/04
// M Fischler	      - put/get to/from streams uses pairs of ulongs when
//			+ storing doubles avoid problems with precision 
//			4/14/05
// Mark Fischler  - Repaired BUG - when mean > 2 billion, was returning
//                  mean instead of the proper value.  01/13/06
// =======================================================================

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/DoubConv.h"
#include <cmath>	// for std::floor()
#include <iostream>
#include <string>
#include <vector>

namespace CLHEP {

std::string RandPoisson::name() const {return "RandPoisson";}
HepRandomEngine & RandPoisson::engine() {return *localEngine;}

// Initialisation of static data
CLHEP_THREAD_LOCAL double RandPoisson::status_st[3] = {0., 0., 0.};
CLHEP_THREAD_LOCAL double RandPoisson::oldm_st = -1.0;
const double RandPoisson::meanMax_st = 2.0E9;

RandPoisson::~RandPoisson() {
}

double RandPoisson::operator()() {
  return double(fire( defaultMean ));
}

double RandPoisson::operator()( double mean ) {
  return double(fire( mean ));
}

double gammln(double xx) {

// Returns the value ln(Gamma(xx) for xx > 0.  Full accuracy is obtained for 
// xx > 1. For 0 < xx < 1. the reflection formula (6.1.4) can be used first.
// (Adapted from Numerical Recipes in C)

  static const double cof[6] = {76.18009172947146,-86.50532032941677,
                             24.01409824083091, -1.231739572450155,
                             0.1208650973866179e-2, -0.5395239384953e-5};
  int j;
  double x = xx - 1.0;
  double tmp = x + 5.5;
  tmp -= (x + 0.5) * std::log(tmp);
  double ser = 1.000000000190015;

  for ( j = 0; j <= 5; j++ ) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp + std::log(2.5066282746310005*ser);
}

static
double normal (HepRandomEngine* eptr) 		// mf 1/13/06
{
  double r;
  double v1,v2,fac;
  do {
    v1 = 2.0 * eptr->flat() - 1.0;
    v2 = 2.0 * eptr->flat() - 1.0;
    r = v1*v1 + v2*v2;
  } while ( r > 1.0 );

  fac = std::sqrt(-2.0*std::log(r)/r);
  return v2*fac;
}

long RandPoisson::shoot(double xm) {

// Returns as a floating-point number an integer value that is a random
// deviation drawn from a Poisson distribution of mean xm, using flat()
// as a source of uniform random numbers.
// (Adapted from Numerical Recipes in C)

  double em, t, y;
  double sq, alxm, g1;
  double om = getOldMean();
  HepRandomEngine* anEngine = HepRandom::getTheEngine();

  double* pstatus = getPStatus();
  sq = pstatus[0];
  alxm = pstatus[1];
  g1 = pstatus[2];

  if( xm == -1 ) return 0;
  if( xm < 12.0 ) {
    if( xm != om ) {
      setOldMean(xm);
      g1 = std::exp(-xm);
    }
    em = -1;
    t = 1.0;
    do {
      em += 1.0;
      t *= anEngine->flat();
    } while( t > g1 );
  }
  else if ( xm < getMaxMean() ) {
    if ( xm != om ) {
      setOldMean(xm);
      sq = std::sqrt(2.0*xm);
      alxm = std::log(xm);
      g1 = xm*alxm - gammln(xm + 1.0);
    }
    do {
      do {
	y = std::tan(CLHEP::pi*anEngine->flat());
	em = sq*y + xm;
      } while( em < 0.0 );
      em = std::floor(em);
      t = 0.9*(1.0 + y*y)* std::exp(em*alxm - gammln(em + 1.0) - g1);
    } while( anEngine->flat() > t );
  }
  else {
    em = xm + std::sqrt(xm) * normal (anEngine);	// mf 1/13/06
    if ( static_cast<long>(em) < 0 ) 
      em = static_cast<long>(xm) >= 0 ? xm : getMaxMean();
  }    
  setPStatus(sq,alxm,g1);
  return long(em);
}

void RandPoisson::shootArray(const int size, long* vect, double m1)
{
  for( long* v = vect; v != vect + size; ++v )
    *v = shoot(m1);
}

long RandPoisson::shoot(HepRandomEngine* anEngine, double xm) {

// Returns as a floating-point number an integer value that is a random
// deviation drawn from a Poisson distribution of mean xm, using flat()
// of a given Random Engine as a source of uniform random numbers.
// (Adapted from Numerical Recipes in C)

  double em, t, y;
  double sq, alxm, g1;
  double om = getOldMean();

  double* pstatus = getPStatus();
  sq = pstatus[0];
  alxm = pstatus[1];
  g1 = pstatus[2];

  if( xm == -1 ) return 0;
  if( xm < 12.0 ) {
    if( xm != om ) {
      setOldMean(xm);
      g1 = std::exp(-xm);
    }
    em = -1;
    t = 1.0;
    do {
      em += 1.0;
      t *= anEngine->flat();
    } while( t > g1 );
  }
  else if ( xm < getMaxMean() ) {
    if ( xm != om ) {
      setOldMean(xm);
      sq = std::sqrt(2.0*xm);
      alxm = std::log(xm);
      g1 = xm*alxm - gammln(xm + 1.0);
    }
    do {
      do {
	y = std::tan(CLHEP::pi*anEngine->flat());
	em = sq*y + xm;
      } while( em < 0.0 );
      em = std::floor(em);
      t = 0.9*(1.0 + y*y)* std::exp(em*alxm - gammln(em + 1.0) - g1);
    } while( anEngine->flat() > t );
  }
  else {
    em = xm + std::sqrt(xm) * normal (anEngine);	// mf 1/13/06
    if ( static_cast<long>(em) < 0 ) 
      em = static_cast<long>(xm) >= 0 ? xm : getMaxMean();
  }    
  setPStatus(sq,alxm,g1);
  return long(em);
}

void RandPoisson::shootArray(HepRandomEngine* anEngine, const int size,
                             long* vect, double m1)
{
  for( long* v = vect; v != vect + size; ++v )
    *v = shoot(anEngine,m1);
}

long RandPoisson::fire() {
  return long(fire( defaultMean ));
}

long RandPoisson::fire(double xm) {

// Returns as a floating-point number an integer value that is a random
// deviation drawn from a Poisson distribution of mean xm, using flat()
// as a source of uniform random numbers.
// (Adapted from Numerical Recipes in C)

  double em, t, y;
  double sq, alxm, g1;

  sq = status[0];
  alxm = status[1];
  g1 = status[2];

  if( xm == -1 ) return 0;
  if( xm < 12.0 ) {
    if( xm != oldm ) {
      oldm = xm;
      g1 = std::exp(-xm);
    }
    em = -1;
    t = 1.0;
    do {
      em += 1.0;
      t *= localEngine->flat();
    } while( t > g1 );
  }
  else if ( xm < meanMax ) {
    if ( xm != oldm ) {
      oldm = xm;
      sq = std::sqrt(2.0*xm);
      alxm = std::log(xm);
      g1 = xm*alxm - gammln(xm + 1.0);
    }
    do {
      do {
	y = std::tan(CLHEP::pi*localEngine->flat());
	em = sq*y + xm;
      } while( em < 0.0 );
      em = std::floor(em);
      t = 0.9*(1.0 + y*y)* std::exp(em*alxm - gammln(em + 1.0) - g1);
    } while( localEngine->flat() > t );
  }
  else {
    em = xm + std::sqrt(xm) * normal (localEngine.get());	// mf 1/13/06
    if ( static_cast<long>(em) < 0 ) 
      em = static_cast<long>(xm) >= 0 ? xm : getMaxMean();
  }    
  status[0] = sq; status[1] = alxm; status[2] = g1;
  return long(em);
}

void RandPoisson::fireArray(const int size, long* vect )
{
  for( long* v = vect; v != vect + size; ++v )
    *v = fire( defaultMean );
}

void RandPoisson::fireArray(const int size, long* vect, double m1)
{
  for( long* v = vect; v != vect + size; ++v )
    *v = fire( m1 );
}

std::ostream & RandPoisson::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  t = DoubConv::dto2longs(meanMax);
  os << meanMax << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(defaultMean);
  os << defaultMean << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(status[0]);
  os << status[0] << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(status[1]);
  os << status[1] << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(status[2]);
  os << status[2] << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(oldm);
  os << oldm << " " << t[0] << " " << t[1] << "\n";
  os.precision(pr);
  return os;
}

std::istream & RandPoisson::get ( std::istream & is ) {
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
  if (possibleKeywordInput(is, "Uvec", meanMax)) {
    std::vector<unsigned long> t(2);
    is >> meanMax     >> t[0] >> t[1]; meanMax     = DoubConv::longs2double(t); 
    is >> defaultMean >> t[0] >> t[1]; defaultMean = DoubConv::longs2double(t); 
    is >> status[0]   >> t[0] >> t[1]; status[0]   = DoubConv::longs2double(t); 
    is >> status[1]   >> t[0] >> t[1]; status[1]   = DoubConv::longs2double(t); 
    is >> status[2]   >> t[0] >> t[1]; status[2]   = DoubConv::longs2double(t); 
    is >> oldm        >> t[0] >> t[1]; oldm        = DoubConv::longs2double(t); 
    return is;
  }
  // is >> meanMax encompassed by possibleKeywordInput
  is >> defaultMean >> status[0] >> status[1] >> status[2];
  return is;
}

}  // namespace CLHEP

