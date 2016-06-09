// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                       --- RandBreitWigner ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
//                - Added methods to shoot arrays: 28th July 1997
// J.Marraffino   - Added default arguments as attributes and
//                  operator() with arguments: 16th Feb 1998
// M Fischler     - put and get to/from streams 12/10/04
// M Fischler	      - put/get to/from streams uses pairs of ulongs when
//			+ storing doubles avoid problems with precision 
//			4/14/05
// =======================================================================

#include "CLHEP/Random/RandBreitWigner.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/DoubConv.h"
#include <algorithm>	// for min() and max()
#include <cmath>

namespace CLHEP {

std::string RandBreitWigner::name() const {return "RandBreitWigner";}
HepRandomEngine & RandBreitWigner::engine() {return *localEngine;}

RandBreitWigner::~RandBreitWigner() {
}

double RandBreitWigner::operator()() {
   return fire( defaultA, defaultB );
}

double RandBreitWigner::operator()( double a, double b ) {
   return fire( a, b );
}

double RandBreitWigner::operator()( double a, double b, double c ) {
   return fire( a, b, c );
}

double RandBreitWigner::shoot(double mean, double gamma)
{
   double rval, displ;

   rval = 2.0*HepRandom::getTheEngine()->flat()-1.0;
   displ = 0.5*gamma*std::tan(rval*CLHEP::halfpi);

   return mean + displ;
}

double RandBreitWigner::shoot(double mean, double gamma, double cut)
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = std::atan(2.0*cut/gamma);
   rval = 2.0*HepRandom::getTheEngine()->flat()-1.0;
   displ = 0.5*gamma*std::tan(rval*val);

   return mean + displ;
}

double RandBreitWigner::shootM2(double mean, double gamma )
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = std::atan(-mean/gamma);
   rval = RandFlat::shoot(val, CLHEP::halfpi);
   displ = gamma*std::tan(rval);

   return std::sqrt(mean*mean + mean*displ);
}

double RandBreitWigner::shootM2(double mean, double gamma, double cut )
{
   double rval, displ;
   double lower, upper, tmp;

   if ( gamma == 0.0 ) return mean;
   tmp = std::max(0.0,(mean-cut));
   lower = std::atan( (tmp*tmp-mean*mean)/(mean*gamma) );
   upper = std::atan( ((mean+cut)*(mean+cut)-mean*mean)/(mean*gamma) );
   rval = RandFlat::shoot(lower, upper);
   displ = gamma*std::tan(rval);

   return std::sqrt(std::max(0.0, mean*mean + mean*displ));
}

void RandBreitWigner::shootArray ( const int size, double* vect )
{
  for( double* v = vect; v != vect + size; ++v )
    *v = shoot( 1.0, 0.2 );
}

void RandBreitWigner::shootArray ( const int size, double* vect,
                                   double a, double b )
{
  for( double* v = vect; v != vect + size; ++v )
    *v = shoot( a, b );
}

void RandBreitWigner::shootArray ( const int size, double* vect,
                                   double a, double b,
                                   double c )
{
  for( double* v = vect; v != vect + size; ++v )
    *v = shoot( a, b, c );
}

//----------------

double RandBreitWigner::shoot(HepRandomEngine* anEngine,
                                 double mean, double gamma)
{
   double rval, displ;

   rval = 2.0*anEngine->flat()-1.0;
   displ = 0.5*gamma*std::tan(rval*CLHEP::halfpi);

   return mean + displ;
}

double RandBreitWigner::shoot(HepRandomEngine* anEngine,
                                 double mean, double gamma, double cut )
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = std::atan(2.0*cut/gamma);
   rval = 2.0*anEngine->flat()-1.0;
   displ = 0.5*gamma*std::tan(rval*val);

   return mean + displ;
}

double RandBreitWigner::shootM2(HepRandomEngine* anEngine,
                                   double mean, double gamma )
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = std::atan(-mean/gamma);
   rval = RandFlat::shoot(anEngine,val, CLHEP::halfpi);
   displ = gamma*std::tan(rval);

   return std::sqrt(mean*mean + mean*displ);
}

double RandBreitWigner::shootM2(HepRandomEngine* anEngine,
                                   double mean, double gamma, double cut )
{
   double rval, displ;
   double lower, upper, tmp;

   if ( gamma == 0.0 ) return mean;
   tmp = std::max(0.0,(mean-cut));
   lower = std::atan( (tmp*tmp-mean*mean)/(mean*gamma) );
   upper = std::atan( ((mean+cut)*(mean+cut)-mean*mean)/(mean*gamma) );
   rval = RandFlat::shoot(anEngine, lower, upper);
   displ = gamma*std::tan(rval);

   return std::sqrt( std::max(0.0, mean*mean+mean*displ) );
}

void RandBreitWigner::shootArray ( HepRandomEngine* anEngine,
                                   const int size, double* vect )
{
  for( double* v = vect; v != vect + size; ++v )
    *v = shoot( anEngine, 1.0, 0.2 );
}

void RandBreitWigner::shootArray ( HepRandomEngine* anEngine,
                                   const int size, double* vect,
                                   double a, double b )
{
  for( double* v = vect; v != vect + size; ++v )
    *v = shoot( anEngine, a, b );
}

void RandBreitWigner::shootArray ( HepRandomEngine* anEngine,
                                   const int size, double* vect,
                                   double a, double b, double c )
{
  for( double* v = vect; v != vect + size; ++v )
    *v = shoot( anEngine, a, b, c );
}

//----------------

double RandBreitWigner::fire()
{
  return fire( defaultA, defaultB );
}

double RandBreitWigner::fire(double mean, double gamma)
{
   double rval, displ;

   rval = 2.0*localEngine->flat()-1.0;
   displ = 0.5*gamma*std::tan(rval*CLHEP::halfpi);

   return mean + displ;
}

double RandBreitWigner::fire(double mean, double gamma, double cut)
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = std::atan(2.0*cut/gamma);
   rval = 2.0*localEngine->flat()-1.0;
   displ = 0.5*gamma*std::tan(rval*val);

   return mean + displ;
}

double RandBreitWigner::fireM2()
{
  return fireM2( defaultA, defaultB );
}

double RandBreitWigner::fireM2(double mean, double gamma )
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = std::atan(-mean/gamma);
   rval = RandFlat::shoot(localEngine.get(),val, CLHEP::halfpi);
   displ = gamma*std::tan(rval);

   return std::sqrt(mean*mean + mean*displ);
}

double RandBreitWigner::fireM2(double mean, double gamma, double cut )
{
   double rval, displ;
   double lower, upper, tmp;

   if ( gamma == 0.0 ) return mean;
   tmp = std::max(0.0,(mean-cut));
   lower = std::atan( (tmp*tmp-mean*mean)/(mean*gamma) );
   upper = std::atan( ((mean+cut)*(mean+cut)-mean*mean)/(mean*gamma) );
   rval = RandFlat::shoot(localEngine.get(),lower, upper);
   displ = gamma*std::tan(rval);

   return std::sqrt(std::max(0.0, mean*mean + mean*displ));
}

void RandBreitWigner::fireArray ( const int size, double* vect)
{
  for( double* v = vect; v != vect + size; ++v )
    *v = fire(defaultA, defaultB );
}

void RandBreitWigner::fireArray ( const int size, double* vect,
                                  double a, double b )
{
  for( double* v = vect; v != vect + size; ++v )
    *v = fire( a, b );
}

void RandBreitWigner::fireArray ( const int size, double* vect,
                                  double a, double b, double c )
{
  for( double* v = vect; v != vect + size; ++v )
    *v = fire( a, b, c );
}


std::ostream & RandBreitWigner::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  t = DoubConv::dto2longs(defaultA);
  os << defaultA << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(defaultB);
  os << defaultB << " " << t[0] << " " << t[1] << "\n";
  os.precision(pr);
  return os;
}

std::istream & RandBreitWigner::get ( std::istream & is ) {
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
    is >> defaultB >> t[0] >> t[1]; defaultB = DoubConv::longs2double(t); 
    return is;
  }
  // is >> defaultA encompassed by possibleKeywordInput
  is >> defaultB;
  return is;
}


}  // namespace CLHEP

