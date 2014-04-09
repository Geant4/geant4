// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                       --- HepJamesRandom ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).
//
// This algorithm implements the original universal random number generator
// as proposed by Marsaglia & Zaman in report FSU-SCRI-87-50 and coded
// in FORTRAN77 by Fred James as the RANMAR generator, part of the MATHLIB
// HEP library.

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
//                - Fixed a bug in setSeed(): 26th February 1996
//                - Minor corrections: 31st October 1996
//                - Added methods for engine status: 19th November 1996
//                - Fixed bug in setSeeds(): 15th September 1997
// J.Marraffino   - Added stream operators and related constructor.
//                  Added automatic seed selection from seed table and
//                  engine counter: 16th Feb 1998
// Ken Smith      - Added conversion operators:  6th Aug 1998
// J. Marraffino  - Remove dependence on hepString class  13 May 1999
// V. Innocente   - changed pointers to indices     3 may 2000
// M. Fischler    - In restore, checkFile for file not found    03 Dec 2004
// M. Fischler    - Methods for distrib. instacne save/restore  12/8/04    
// M. Fischler    - split get() into tag validation and 
//                  getState() for anonymous restores           12/27/04    
// M. Fischler    - Enforcement that seeds be non-negative
//		    (lest the sequence be non-random)	         2/14/05    
// M. Fischler    - put/get for vectors of ulongs		3/14/05
// M. Fischler    - State-saving using only ints, for portability 4/12/05
//		    
// =======================================================================

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/engineIDulong.h"
#include "CLHEP/Random/DoubConv.h"
#include "CLHEP/Utility/atomic_int.h"

#include <string.h>	// for strcmp
#include <cmath>
#include <cstdlib>

namespace CLHEP {

namespace {
  // Number of instances with automatic seed selection
  CLHEP_ATOMIC_INT_TYPE numberOfEngines(0);

  // Maximum index into the seed table
  const int maxIndex = 215;
}

static const int MarkerLen = 64; // Enough room to hold a begin or end marker. 

std::string HepJamesRandom::name() const {return "HepJamesRandom";}

HepJamesRandom::HepJamesRandom(long seed)
: HepRandomEngine()
{
  setSeed(seed,0);
  setSeeds(&theSeed,0);
}

HepJamesRandom::HepJamesRandom()     	// 15 Feb. 1998  JMM
: HepRandomEngine()
{
  long seeds[2];
  long seed;

  int numEngines = numberOfEngines++;
  int cycle = std::abs(int(numEngines/maxIndex));
  int curIndex = std::abs(int(numEngines%maxIndex));

  long mask = ((cycle & 0x007fffff) << 8);
  HepRandom::getTheTableSeeds( seeds, curIndex );
  seed = seeds[0]^mask;
  setSeed(seed,0);
  setSeeds(&theSeed,0);
}

HepJamesRandom::HepJamesRandom(int rowIndex, int colIndex) // 15 Feb. 1998  JMM
: HepRandomEngine()
{
  long seed;
   long seeds[2];

  int cycle = std::abs(int(rowIndex/maxIndex));
  int row = std::abs(int(rowIndex%maxIndex));
  int col = std::abs(int(colIndex%2));
  long mask = ((cycle & 0x000007ff) << 20);
  HepRandom::getTheTableSeeds( seeds, row );
  seed = (seeds[col])^mask;
  setSeed(seed,0);
  setSeeds(&theSeed,0);
}

HepJamesRandom::HepJamesRandom(std::istream& is)
: HepRandomEngine()
{
  is >> *this;
}

HepJamesRandom::~HepJamesRandom() {}

void HepJamesRandom::saveStatus( const char filename[] ) const
{
  std::ofstream outFile( filename, std::ios::out ) ;

  if (!outFile.bad()) {
    outFile << "Uvec\n";
    std::vector<unsigned long> v = put();
    for (unsigned int i=0; i<v.size(); ++i) {
      outFile << v[i] << "\n";
    }
  }
}

void HepJamesRandom::restoreStatus( const char filename[] )
{
   int ipos, jpos;
   std::ifstream inFile( filename, std::ios::in);
   if (!checkFile ( inFile, filename, engineName(), "restoreStatus" )) {
     std::cerr << "  -- Engine state remains unchanged\n";
     return;
   }
  if ( possibleKeywordInput ( inFile, "Uvec", theSeed ) ) {
    std::vector<unsigned long> v;
    unsigned long xin;
    for (unsigned int ivec=0; ivec < VECTOR_STATE_SIZE; ++ivec) {
      inFile >> xin;
      if (!inFile) {
        inFile.clear(std::ios::badbit | inFile.rdstate());
        std::cerr << "\nJamesRandom state (vector) description improper."
	       << "\nrestoreStatus has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return;
      }
      v.push_back(xin);
    }
    getState(v);
    return;
  }

   if (!inFile.bad() && !inFile.eof()) {
//     inFile >> theSeed;  removed -- encompased by possibleKeywordInput
     for (int i=0; i<97; ++i)
       inFile >> u[i];
     inFile >> c; inFile >> cd; inFile >> cm;
     inFile >> jpos;
     ipos = (64+jpos)%97;
     i97 = ipos;
     j97 = jpos;
   }
}

void HepJamesRandom::showStatus() const
{
   std::cout << std::endl;
   std::cout << "----- HepJamesRandom engine status -----" << std::endl;
   std::cout << " Initial seed = " << theSeed << std::endl;
   std::cout << " u[] = ";
   for (int i=0; i<97; ++i)
     std::cout << u[i] << " ";
   std::cout << std::endl;
   std::cout << " c = " << c << ", cd = " << cd << ", cm = " << cm
	     << std::endl;
   std::cout << " i97 = " << i97 << ", u[i97] = " << u[i97] << std::endl;
   std::cout << " j97 = " << j97 << ", u[j97] = " << u[j97] << std::endl;
   std::cout << "----------------------------------------" << std::endl;
}

void HepJamesRandom::setSeed(long seed, int)
{
  // The input value for "seed" should be within the range [0,900000000]
  //
  // Negative seeds result in serious flaws in the randomness;
  // seeds above 900000000 are OK because of the %177 in the expression for i,
  // but may have the same effect as other seeds below 900000000.

  int m, n;
  float s, t;
  long mm;

  if (seed < 0) {
    std::cout << "Seed for HepJamesRandom must be non-negative\n" 
    	<< "Seed value supplied was " << seed  
	<< "\nUsing its absolute value instead\n";
    seed = -seed;
  }
  
  long ij = seed/30082;
  long kl = seed - 30082*ij;
  long i = (ij/177) % 177 + 2;
  long j = ij % 177 + 2;
  long k = (kl/169) % 178 + 1;
  long l = kl % 169;

  theSeed = seed;

  for ( n = 1 ; n < 98 ; n++ ) {
    s = 0.0;
    t = 0.5;
    for ( m = 1 ; m < 25 ; m++) {
      mm = ( ( (i*j) % 179 ) * k ) % 179;
      i = j;
      j = k;
      k = mm;
      l = ( 53 * l + 1 ) % 169;
      if ( (l*mm % 64 ) >= 32 )
        s += t;
      t *= 0.5;
    }
    u[n-1] = s;
  }
  c = 362436.0 / 16777216.0;
  cd = 7654321.0 / 16777216.0;
  cm = 16777213.0 / 16777216.0;

  i97 = 96;
  j97 = 32;

}

void HepJamesRandom::setSeeds(const long* seeds, int)
{
  setSeed(seeds ? *seeds : 19780503L, 0);
  theSeeds = seeds;
}

double HepJamesRandom::flat()
{
   double uni;

   do {
      uni = u[i97] - u[j97];
      if ( uni < 0.0 ) uni++;
      u[i97] = uni;
      
      if (i97 == 0) i97 = 96;
      else i97--;
      
      if (j97 == 0) j97 = 96;
      else j97--;
      
      c -= cd;
      if (c < 0.0) c += cm;
      
      uni -= c;
      if (uni < 0.0) uni += 1.0;
   } while ( uni <= 0.0 || uni >= 1.0 );
   
   return uni;
}

void HepJamesRandom::flatArray(const int size, double* vect)
{
//   double uni;
   int i;

   for (i=0; i<size; ++i) {
     vect[i] = flat();
   }   
}

HepJamesRandom::operator unsigned int() {
   return ((unsigned int)(flat() * exponent_bit_32()) & 0xffffffff )  |
         (((unsigned int)( u[i97] * exponent_bit_32())>>16)  & 0xff);
}

std::ostream & HepJamesRandom::put ( std::ostream& os ) const {
  char beginMarker[] = "JamesRandom-begin";
  os << beginMarker << "\nUvec\n";
  std::vector<unsigned long> v = put();
  for (unsigned int i=0; i<v.size(); ++i) {
     os <<  v[i] <<  "\n";
  }
  return os;  
}

std::vector<unsigned long> HepJamesRandom::put () const {
  std::vector<unsigned long> v;
  v.push_back (engineIDulong<HepJamesRandom>());
  std::vector<unsigned long> t;
  for (int i=0; i<97; ++i) {
    t = DoubConv::dto2longs(u[i]);
    v.push_back(t[0]); v.push_back(t[1]);
  }
  t = DoubConv::dto2longs(c);
  v.push_back(t[0]); v.push_back(t[1]);
  t = DoubConv::dto2longs(cd);
  v.push_back(t[0]); v.push_back(t[1]);
  t = DoubConv::dto2longs(cm);
  v.push_back(t[0]); v.push_back(t[1]);
  v.push_back(static_cast<unsigned long>(j97));
  return v;
}


std::istream & HepJamesRandom::get  ( std::istream& is) {
  char beginMarker [MarkerLen];
  is >> std::ws;
  is.width(MarkerLen);  // causes the next read to the char* to be <=
			// that many bytes, INCLUDING A TERMINATION \0 
			// (Stroustrup, section 21.3.2)
  is >> beginMarker;
  if (strcmp(beginMarker,"JamesRandom-begin")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nInput stream mispositioned or"
	       << "\nJamesRandom state description missing or"
	       << "\nwrong engine type found." << std::endl;
     return is;
  }
  return getState(is);
}

std::string HepJamesRandom::beginTag ( )  { 
  return "JamesRandom-begin"; 
}

std::istream & HepJamesRandom::getState  ( std::istream& is) {
  if ( possibleKeywordInput ( is, "Uvec", theSeed ) ) {
    std::vector<unsigned long> v;
    unsigned long uu;
    for (unsigned int ivec=0; ivec < VECTOR_STATE_SIZE; ++ivec) {
      is >> uu;
      if (!is) {
        is.clear(std::ios::badbit | is.rdstate());
        std::cerr << "\nJamesRandom state (vector) description improper."
		<< "\ngetState() has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return is;
      }
      v.push_back(uu);
    }
    getState(v);
    return (is);
  }

//  is >> theSeed;  Removed, encompassed by possibleKeywordInput()

  int ipos, jpos;
  char   endMarker [MarkerLen];
  for (int i=0; i<97; ++i) {
     is >> u[i];
  }
  is >> c; is >> cd; is >> cm;
  is >> jpos;
  is >> std::ws;
  is.width(MarkerLen);
  is >> endMarker;
  if(strcmp(endMarker,"JamesRandom-end")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nJamesRandom state description incomplete."
	       << "\nInput stream is probably mispositioned now." << std::endl;
     return is;
  }

  ipos = (64+jpos)%97;
  i97 = ipos;
  j97 = jpos;
  return is;
}

bool HepJamesRandom::get (const std::vector<unsigned long> & v) {
  if ( (v[0] & 0xffffffffUL) != engineIDulong<HepJamesRandom>()) {
    std::cerr << 
    	"\nHepJamesRandom get:state vector has wrong ID word - state unchanged\n";
    return false;
  }
  return getState(v);
}

bool HepJamesRandom::getState (const std::vector<unsigned long> & v) {
  if (v.size() != VECTOR_STATE_SIZE ) {
    std::cerr << 
    	"\nHepJamesRandom get:state vector has wrong length - state unchanged\n";
    return false;
  }
  std::vector<unsigned long> t(2);
  for (int i=0; i<97; ++i) {
    t[0] = v[2*i+1]; t[1] = v[2*i+2];
    u[i] = DoubConv::longs2double(t);
  }
  t[0] = v[195]; t[1] = v[196]; c  = DoubConv::longs2double(t);
  t[0] = v[197]; t[1] = v[198]; cd = DoubConv::longs2double(t);
  t[0] = v[199]; t[1] = v[200]; cm = DoubConv::longs2double(t);
  j97  = v[201];
  i97  = (64+j97)%97; 
  return true;
}

}  // namespace CLHEP
