// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                        --- MTwistEngine ---
//                      class implementation file
// -----------------------------------------------------------------------
// A "fast, compact, huge-period generator" based on M. Matsumoto and 
// T. Nishimura, "Mersenne Twister: A 623-dimensionally equidistributed 
// uniform pseudorandom number generator", to appear in ACM Trans. on
// Modeling and Computer Simulation.  It is a twisted GFSR generator
// with a Mersenne-prime period of 2^19937-1, uniform on open interval (0,1)
// =======================================================================
// Ken Smith      - Started initial draft:                    14th Jul 1998
//                - Optimized to get std::pow() out of flat() method
//                - Added conversion operators:                6th Aug 1998
// J. Marraffino  - Added some explicit casts to deal with
//                  machines where sizeof(int) != sizeof(long)  22 Aug 1998
// M. Fischler    - Modified constructors such that no two
//		    seeds will match sequences, no single/double
//		    seeds will match, explicit seeds give 
//		    determined results, and default will not 
//		    match any of the above or other defaults.
//                - Modified use of the various exponents of 2
//                  to avoid per-instance space overhead and
//                  correct the rounding procedure              16 Sep 1998
// J. Marfaffino  - Remove dependence on hepString class        13 May 1999
// M. Fischler    - In restore, checkFile for file not found    03 Dec 2004
// M. Fischler    - Methods for distrib. instacne save/restore  12/8/04    
// M. Fischler    - split get() into tag validation and 
//                  getState() for anonymous restores           12/27/04    
// M. Fischler    - put/get for vectors of ulongs		3/14/05
// M. Fischler    - State-saving using only ints, for portability 4/12/05
// M. Fischler    - Improved seeding in setSeed  (Savanah bug #17479) 11/15/06
//		  - (Possible optimization - now that the seeding is improved,
//		    is it still necessary for ctor to "warm up" by discarding
//		    2000 iterations?)
//		    
// =======================================================================

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/engineIDulong.h"
#include "CLHEP/Utility/atomic_int.h"

#include <string.h>	// for strcmp
#include <cstdlib>	// for std::abs(int)

namespace CLHEP {

namespace {
  // Number of instances with automatic seed selection
  CLHEP_ATOMIC_INT_TYPE numberOfEngines(0);

  // Maximum index into the seed table
  const int maxIndex = 215;
}

static const int MarkerLen = 64; // Enough room to hold a begin or end marker. 

std::string MTwistEngine::name() const {return "MTwistEngine";}

MTwistEngine::MTwistEngine() 
: HepRandomEngine()
{
  int numEngines = numberOfEngines++;
  int cycle = std::abs(int(numEngines/maxIndex));
  int curIndex = std::abs(int(numEngines%maxIndex));
  long mask = ((cycle & 0x007fffff) << 8);
  long seedlist[2];
  HepRandom::getTheTableSeeds( seedlist, curIndex );
  seedlist[0] = (seedlist[0])^mask;
  seedlist[1] = 0;
  setSeeds( seedlist, numEngines );
  count624=0;

  for( int i=0; i < 2000; ++i ) flat();      // Warm up just a bit
}

MTwistEngine::MTwistEngine(long seed)  
: HepRandomEngine()
{
  long seedlist[2]={seed,17587};
  setSeeds( seedlist, 0 );
  count624=0;
  for( int i=0; i < 2000; ++i ) flat();      // Warm up just a bit
}

MTwistEngine::MTwistEngine(int rowIndex, int colIndex) 
: HepRandomEngine()
{
  int cycle = std::abs(int(rowIndex/maxIndex));
  int row = std::abs(int(rowIndex%maxIndex));
  int col = std::abs(int(colIndex%2));
  long mask = (( cycle & 0x000007ff ) << 20 );
  long seedlist[2];
  HepRandom::getTheTableSeeds( seedlist, row );
  seedlist[0] = (seedlist[col])^mask;
  seedlist[1] = 690691;
  setSeeds(seedlist, 4444772);
  count624=0;
  for( int i=0; i < 2000; ++i ) flat();      // Warm up just a bit
}

MTwistEngine::MTwistEngine( std::istream& is )  
: HepRandomEngine()
{
  is >> *this;
}

MTwistEngine::~MTwistEngine() {}

double MTwistEngine::flat() {
  unsigned int y;

   if( count624 >= N ) {
    int i;

    for( i=0; i < NminusM; ++i ) {
      y = (mt[i] & 0x80000000) | (mt[i+1] & 0x7fffffff);
      mt[i] = mt[i+M]       ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );
    }

    for(    ; i < N-1    ; ++i ) {
      y = (mt[i] & 0x80000000) | (mt[i+1] & 0x7fffffff);
      mt[i] = mt[i-NminusM] ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );
    }

    y = (mt[i] & 0x80000000) | (mt[0] & 0x7fffffff);
    mt[i] = mt[M-1] ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );

    count624 = 0;
  }

  y = mt[count624];
  y ^= ( y >> 11);
  y ^= ((y << 7 ) & 0x9d2c5680);
  y ^= ((y << 15) & 0xefc60000);
  y ^= ( y >> 18);

  return                      y * twoToMinus_32()  +    // Scale to range 
         (mt[count624++] >> 11) * twoToMinus_53()  +    // fill remaining bits
                	    nearlyTwoToMinus_54();      // make sure non-zero
}

void MTwistEngine::flatArray( const int size, double *vect ) {
  for( int i=0; i < size; ++i) vect[i] = flat();
}

void MTwistEngine::setSeed(long seed, int k) {

  // MF 11/15/06 - Change seeding algorithm to a newer one recommended 
  //               by Matsumoto: The original algorithm was 
  //		   mt[i] = (69069 * mt[i-1]) & 0xffffffff and this gives
  //		   problems when the seed bit pattern has lots of zeros
  //		   (for example, 0x08000000).  Savanah bug #17479.

  theSeed = seed ? seed : 4357;
  int mti;
  const int N1=624;
  mt[0] = (unsigned int) (theSeed&0xffffffffUL);
  for (mti=1; mti<N1; mti++) {
    mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
    mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
  }
  for( int i=1; i < 624; ++i ) {
    mt[i] ^= k;			// MF 9/16/98: distinguish starting points
  }
  // MF 11/15/06 This distinction of starting points based on values of k
  //             is kept even though the seeding algorithm itself is improved.
}

void MTwistEngine::setSeeds(const long *seeds, int k) {
  setSeed( (*seeds ? *seeds : 43571346), k );
  int i;
  for( i=1; i < 624; ++i ) {
    mt[i] = ( seeds[1] + mt[i] ) & 0xffffffff; // MF 9/16/98 
  }
  theSeeds = seeds;
}

void MTwistEngine::saveStatus( const char filename[] ) const
{
   std::ofstream outFile( filename, std::ios::out ) ;
   if (!outFile.bad()) {
     outFile << theSeed << std::endl;
     for (int i=0; i<624; ++i) outFile <<std::setprecision(20) << mt[i] << " ";
     outFile << std::endl;
     outFile << count624 << std::endl;
   }
}

void MTwistEngine::restoreStatus( const char filename[] )
{
   std::ifstream inFile( filename, std::ios::in);
   if (!checkFile ( inFile, filename, engineName(), "restoreStatus" )) {
     std::cerr << "  -- Engine state remains unchanged\n";
     return;
   }

   if (!inFile.bad() && !inFile.eof()) {
     inFile >> theSeed;
     for (int i=0; i<624; ++i) inFile >> mt[i];
     inFile >> count624;
   }
}

void MTwistEngine::showStatus() const
{
   std::cout << std::endl;
   std::cout << "--------- MTwist engine status ---------" << std::endl;
   std::cout << std::setprecision(20);
   std::cout << " Initial seed      = " << theSeed << std::endl;
   std::cout << " Current index     = " << count624 << std::endl;
   std::cout << " Array status mt[] = " << std::endl;
   // 2014/06/06  L Garren
   // the final line has 4 elements, not 5
   for (int i=0; i<620; i+=5) {
     std::cout << mt[i]   << " " << mt[i+1] << " " << mt[i+2] << " " 
	       << mt[i+3] << " " << mt[i+4] << "\n";
   }
   std::cout << mt[620]   << " " << mt[621] << " " << mt[622] << " " 
	     << mt[623]  << std::endl;
   std::cout << "----------------------------------------" << std::endl;
}

MTwistEngine::operator float() {
  unsigned int y;

  if( count624 >= N ) {
    int i;

    for( i=0; i < NminusM; ++i ) {
      y = (mt[i] & 0x80000000) | (mt[i+1] & 0x7fffffff);
      mt[i] = mt[i+M]       ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );
    }

    for(    ; i < N-1    ; ++i ) {
      y = (mt[i] & 0x80000000) | (mt[i+1] & 0x7fffffff);
      mt[i] = mt[i-NminusM] ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );
    }

    y = (mt[i] & 0x80000000) | (mt[0] & 0x7fffffff);
    mt[i] = mt[M-1] ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );

    count624 = 0;
  }

  y = mt[count624++];
  y ^= ( y >> 11);
  y ^= ((y << 7 ) & 0x9d2c5680);
  y ^= ((y << 15) & 0xefc60000);
  y ^= ( y >> 18);

  return (float)(y * twoToMinus_32());
}

MTwistEngine::operator unsigned int() {
  unsigned int y;

  if( count624 >= N ) {
    int i;

    for( i=0; i < NminusM; ++i ) {
      y = (mt[i] & 0x80000000) | (mt[i+1] & 0x7fffffff);
      mt[i] = mt[i+M]       ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );
    }

    for(    ; i < N-1    ; ++i ) {
      y = (mt[i] & 0x80000000) | (mt[i+1] & 0x7fffffff);
      mt[i] = mt[i-NminusM] ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );
    }

    y = (mt[i] & 0x80000000) | (mt[0] & 0x7fffffff);
    mt[i] = mt[M-1] ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );

    count624 = 0;
  }

  y = mt[count624++];
  y ^= ( y >> 11);
  y ^= ((y << 7 ) & 0x9d2c5680);
  y ^= ((y << 15) & 0xefc60000);
  y ^= ( y >> 18);

  return y;
}

std::ostream & MTwistEngine::put ( std::ostream& os ) const
{
   char beginMarker[] = "MTwistEngine-begin";
   char endMarker[]   = "MTwistEngine-end";

   int pr = os.precision(20);
   os << " " << beginMarker << " ";
   os << theSeed << " ";
   for (int i=0; i<624; ++i) {
     os << mt[i] << "\n";
   }
   os << count624 << " ";
   os << endMarker << "\n";
   os.precision(pr);
   return os;
}

std::vector<unsigned long> MTwistEngine::put () const {
  std::vector<unsigned long> v;
  v.push_back (engineIDulong<MTwistEngine>());
  for (int i=0; i<624; ++i) {
     v.push_back(static_cast<unsigned long>(mt[i]));
  }
  v.push_back(count624); 
  return v;
}

std::istream &  MTwistEngine::get ( std::istream& is )
{
  char beginMarker [MarkerLen];
  is >> std::ws;
  is.width(MarkerLen);  // causes the next read to the char* to be <=
			// that many bytes, INCLUDING A TERMINATION \0 
			// (Stroustrup, section 21.3.2)
  is >> beginMarker;
  if (strcmp(beginMarker,"MTwistEngine-begin")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nInput stream mispositioned or"
	       << "\nMTwistEngine state description missing or"
	       << "\nwrong engine type found." << std::endl;
     return is;
   }
  return getState(is);
}

std::string MTwistEngine::beginTag ( )  { 
  return "MTwistEngine-begin"; 
}

std::istream &  MTwistEngine::getState ( std::istream& is )
{
  char endMarker   [MarkerLen];
  is >> theSeed;
  for (int i=0; i<624; ++i)  is >> mt[i];
  is >> count624;
  is >> std::ws;
  is.width(MarkerLen);  
  is >> endMarker;
  if (strcmp(endMarker,"MTwistEngine-end")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nMTwistEngine state description incomplete."
	       << "\nInput stream is probably mispositioned now." << std::endl;
     return is;
   }
   return is;
}

bool MTwistEngine::get (const std::vector<unsigned long> & v) {
  if ((v[0] & 0xffffffffUL) != engineIDulong<MTwistEngine>()) {
    std::cerr << 
    	"\nMTwistEngine get:state vector has wrong ID word - state unchanged\n";
    return false;
  }
  return getState(v);
}

bool MTwistEngine::getState (const std::vector<unsigned long> & v) {
  if (v.size() != VECTOR_STATE_SIZE ) {
    std::cerr << 
    	"\nMTwistEngine get:state vector has wrong length - state unchanged\n";
    return false;
  }
  for (int i=0; i<624; ++i) {
     mt[i]=v[i+1];
  }
  count624 = v[625];
  return true;
}

}  // namespace CLHEP
