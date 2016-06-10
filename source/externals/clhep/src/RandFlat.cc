// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- RandFlat ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// =======================================================================
// Gabriele Cosmo - Created: 17th May 1995
//                - Added methods to shoot arrays: 28th July 1997
//                - Added operator(): 24th Jul 1997
// J.Marraffino   - Added default arguments as attributes and
//                  operator() with arguments: 16th Feb 1998
// M Fischler     - Copy constructor should supply right engine to HepRandom:
//                  1/26/00.
// M Fischler	  - Semi-fix to the saveEngineStatus misbehavior causing
//		    non-reproducing shootBit() 3/1/00.
// M Fischler     - Avoiding hang when file not found in restoreEngineStatus 
//                  12/3/04
// M Fischler     - put and get to/from streams 12/10/04
// M Fischler     - save and restore dist to streams 12/20/04
// M Fischler	      - put/get to/from streams uses pairs of ulongs when
//			+ storing doubles avoid problems with precision 
//			4/14/05
// =======================================================================

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/DoubConv.h"
#include <string.h>	// for strcmp

namespace CLHEP {

const int RandFlat::MSBBits= 15;
const unsigned long RandFlat::MSB= 1ul<<RandFlat::MSBBits;
CLHEP_THREAD_LOCAL unsigned long RandFlat::staticRandomInt= 0;
CLHEP_THREAD_LOCAL unsigned long RandFlat::staticFirstUnusedBit= 0;

std::string RandFlat::name() const {return "RandFlat";}
HepRandomEngine & RandFlat::engine() {return *localEngine;}

RandFlat::~RandFlat() {
}

double RandFlat::operator()() {
  return fire( defaultA, defaultB );
}

double RandFlat::operator()( double w ) {
  return fire( w );
}

double RandFlat::operator()( double a, double b ) {
  return fire( a, b );
}

double RandFlat::shoot() {
  return HepRandom::getTheEngine()->flat();
}

void RandFlat::shootArray(const int size, double* vect) {
  HepRandom::getTheEngine()->flatArray(size,vect);
}

void RandFlat::shootArray( const int size, double* vect,
                           double lx, double dx  )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot(lx,dx);
}

void RandFlat::shootArray( HepRandomEngine* anEngine,
                           const int size, double* vect,
                           double lx, double dx  )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot(anEngine,lx,dx);
}

void RandFlat::fireArray( const int size, double* vect)
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = fire( defaultA, defaultB );
}

void RandFlat::fireArray( const int size, double* vect,
                          double lx, double dx  )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = fire( lx, dx );
}

void RandFlat::saveEngineStatus ( const char filename[] ) {

  // First save the engine status just like the base class would do:
  getTheEngine()->saveStatus( filename );

  // Now append the cached random Int, and first unused bit:

  std::ofstream outfile ( filename, std::ios::app );

  outfile << "RANDFLAT staticRandomInt: " << staticRandomInt 
          << "    staticFirstUnusedBit: " << staticFirstUnusedBit << "\n";

} // saveEngineStatus


void RandFlat::restoreEngineStatus( const char filename[] ) {

  // First restore the engine status just like the base class would do:
  getTheEngine()->restoreStatus( filename );

  // Now find the line describing the cached data:

  std::ifstream infile ( filename, std::ios::in );
  if (!infile) return;
  char inputword[] = "NO_KEYWORD    "; // leaves room for 14 characters plus \0
  while (true) {
    infile.width(13);
    infile >> inputword;
    if (strcmp(inputword,"RANDFLAT")==0) break;
    if (infile.eof()) break;
        // If the file ends without the RANDFLAT line, that means this
        // was a file produced by an earlier version of RandFlat.  We will
        // replicate the old behavior in that case:  staticFirstUnusedBit 
	// and staticRandomInt retain their existing values.
  }

  // Then read and use the caching info:

  if (strcmp(inputword,"RANDFLAT")==0) {
    char setword[40];	// the longest, staticFirstUnusedBit: has length 21
    infile.width(39);
    infile >> setword;
    // setword should be staticRandomInt:  
    infile >> staticRandomInt;
    infile.width(39);
    infile >> setword;
    // setword should be staticFirstUnusedBit: 
    infile >> staticFirstUnusedBit;
  }

} // restoreEngineStatus

std::ostream & RandFlat::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  os << randomInt << " " << firstUnusedBit << "\n";
  t = DoubConv::dto2longs(defaultWidth);
  os << defaultWidth << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(defaultA);
  os << defaultA << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(defaultB);
  os << defaultB << " " << t[0] << " " << t[1] << "\n";
  os.precision(pr);
  return os;
}

std::istream & RandFlat::get ( std::istream & is ) {
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
  if (possibleKeywordInput(is, "Uvec", randomInt)) {
    std::vector<unsigned long> t(2);
    is >> randomInt >> firstUnusedBit;
    is >> defaultWidth >>t[0]>>t[1]; defaultWidth = DoubConv::longs2double(t); 
    is >> defaultA >> t[0] >> t[1]; defaultA = DoubConv::longs2double(t); 
    is >> defaultB >> t[0] >> t[1]; defaultB = DoubConv::longs2double(t); 
    if (!is) {
      is.clear(std::ios::badbit | is.rdstate());
      std::cerr << "\nRandFlat input failed"
	     << "\nInput stream is probably mispositioned now." << std::endl;
      return is;
    }
    return is;
  }
  // is >> randomInt encompassed by possibleKeywordInput
  is >> firstUnusedBit;
  is >> defaultWidth >> defaultA >> defaultB;
  return is;
}

std::ostream & RandFlat::saveDistState ( std::ostream & os ) {
  os << distributionName() << "\n";
  int prec = os.precision(20);
  os << "RANDFLAT staticRandomInt: " << staticRandomInt 
     << "    staticFirstUnusedBit: " << staticFirstUnusedBit << "\n";
  os.precision(prec);
  return os;
}    

std::istream & RandFlat::restoreDistState ( std::istream & is ) {
  std::string inName;
  is >> inName;
  if (inName != distributionName()) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Mismatch when expecting to read static state of a "
    	      << distributionName() << " distribution\n"
	      << "Name found was " << inName
	      << "\nistream is left in the badbit state\n";
    return is;
  }
  std::string keyword;
  std::string c1;
  std::string c2;
  is >> keyword;
  if (keyword!="RANDFLAT") {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Mismatch when expecting to read RANDFLAT bit cache info: "
    	      << keyword << "\n";
    return is;
  }
  is >> c1 >> staticRandomInt >> c2 >> staticFirstUnusedBit;
  return is;
} 

std::ostream & RandFlat::saveFullState ( std::ostream & os ) {
  HepRandom::saveFullState(os);
  saveDistState(os);
  return os;
}
  
std::istream & RandFlat::restoreFullState ( std::istream & is ) {
  HepRandom::restoreFullState(is);
  restoreDistState(is);
  return is;
}


}  // namespace CLHEP

