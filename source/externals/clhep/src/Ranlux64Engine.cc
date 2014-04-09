// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                       --- Ranlux64Engine ---
//                      class implementation file
// -----------------------------------------------------------------------
// A double-precision implementation of the RanluxEngine generator as 
// decsribed by the notes of the original ranlux author (Martin Luscher)
//
// See the note by Martin Luscher, December 1997, entitiled
// Double-precision implementation of the random number generator ranlux
//
// =======================================================================
// Ken Smith      - Initial draft: 14th Jul 1998
//                - Removed pow() from flat method 14th Jul 1998
//                - Added conversion operators:  6th Aug 1998
//
// Mark Fischler  The following were modified mostly to make the routine
//		  exactly match the Luscher algorithm in generating 48-bit
//		  randoms:
// 9/9/98	  - Substantial changes in what used to be flat() to match
//		    algorithm in Luscher's ranlxd.c
//		  - Added update() method for 12 numbers, making flat() trivial
//		  - Added advance() method to hold the unrolled loop for update
//		  - Distinction between three forms of seeding such that it
//		    is impossible to get same sequence from different forms -
//		    done by discarding some fraction of one macro cycle which
//		    is different for the three cases
//		  - Change the misnomer "seed_table" to the more accurate 
//		    "randoms"
//		  - Removed the no longer needed count12, i_lag, j_lag, etc.
//		  - Corrected seed procedure which had been filling bits past
//		    2^-48.  This actually was very bad, invalidating the
//		    number theory behind the proof that ranlxd is good.
//		  - Addition of 2**(-49) to generated number to prevent zero 
//		    from being returned; this does not affect the sequence 
//		    itself.
//		  - Corrected ecu seeding, which had been supplying only 
//		    numbers less than 1/2.  This is probably moot.
// 9/15/98	  - Modified use of the various exponents of 2
//                  to avoid per-instance space overhead.  Note that these
//		    are initialized in setSeed, which EVERY constructor
//		    must invoke.
// J. Marraffino  - Remove dependence on hepString class  13 May 1999
// M. Fischler    - In restore, checkFile for file not found    03 Dec 2004
// M. Fischler    - put get Methods for distrib instance save/restore 12/8/04    
// M. Fischler    - split get() into tag validation and 
//                  getState() for anonymous restores           12/27/04    
// M. Fischler    - put/get for vectors of ulongs		3/14/05
// M. Fischler    - State-saving using only ints, for portability 4/12/05
//
// =======================================================================

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include "CLHEP/Random/engineIDulong.h"
#include "CLHEP/Random/DoubConv.h"
#include "CLHEP/Utility/atomic_int.h"

#include <string.h>	// for strcmp
#include <cstdlib>	// for std::abs(int)
#include <limits>       // for numeric_limits

namespace CLHEP {

namespace {
  // Number of instances with automatic seed selection
  CLHEP_ATOMIC_INT_TYPE numberOfEngines(0);

  // Maximum index into the seed table
  const int maxIndex = 215;
}

static const int MarkerLen = 64; // Enough room to hold a begin or end marker. 


#ifndef WIN32
namespace detail {

template< std::size_t n,
          bool = n < std::size_t(std::numeric_limits<unsigned long>::digits) >
  struct do_right_shift;
template< std::size_t n >
  struct do_right_shift<n,true>
{
  unsigned long operator()(unsigned long value) { return value >> n; }
};
template< std::size_t n >
  struct do_right_shift<n,false>
{
  unsigned long operator()(unsigned long) { return 0ul; }
};

template< std::size_t nbits >
  unsigned long rshift( unsigned long value )
{ return do_right_shift<nbits>()(value); }

} // namespace detail
#endif

std::string Ranlux64Engine::name() const {return "Ranlux64Engine";}

Ranlux64Engine::Ranlux64Engine()
: HepRandomEngine()
{
   luxury = 1;
   int numEngines = numberOfEngines++;
   int cycle    = std::abs(int(numEngines/maxIndex));
   int curIndex = std::abs(int(numEngines%maxIndex));

   long mask = ((cycle & 0x007fffff) << 8);
   long seedlist[2];
   HepRandom::getTheTableSeeds( seedlist, curIndex );
   seedlist[0] ^= mask;
   seedlist[1] = 0;

   setSeeds(seedlist, luxury);
   advance ( 8 );  		// Discard some iterations and ensure that
				// this sequence won't match one where seeds 
				// were provided.
}

Ranlux64Engine::Ranlux64Engine(long seed, int lux)
: HepRandomEngine()
{
   luxury = lux;
   long seedlist[2]={seed,0};
   setSeeds(seedlist, lux);
   advance ( 2*lux + 1 );  	// Discard some iterations to use a different 
				// point in the sequence.  
}

Ranlux64Engine::Ranlux64Engine(int rowIndex, int, int lux)
: HepRandomEngine()
{
   luxury = lux;
   int cycle = std::abs(int(rowIndex/maxIndex));
   int   row = std::abs(int(rowIndex%maxIndex));
   long mask = (( cycle & 0x000007ff ) << 20 );
   long seedlist[2]; 
   HepRandom::getTheTableSeeds( seedlist, row );
   seedlist[0] ^= mask;
   seedlist[1]= 0;
   setSeeds(seedlist, lux);
}

Ranlux64Engine::Ranlux64Engine( std::istream& is )
: HepRandomEngine()
{
  is >> *this;
}

Ranlux64Engine::~Ranlux64Engine() {}

double Ranlux64Engine::flat() {
  // Luscher improves the speed by computing several numbers in a shot,
  // in a manner similar to that of the Tausworth in DualRand or the Hurd
  // engines.  Thus, the real work is done in update().  Here we merely ensure
  // that zero, which the algorithm can produce, is never returned by flat().

  if (index <= 0) update();
  return randoms[--index] + twoToMinus_49();
}

void Ranlux64Engine::update() {
  // Update the stash of twelve random numbers.  
  // When this routione is entered, index is always 0.  The randoms 
  // contains the last 12 numbers in the sequents:  s[0] is x[a+11], 
  // s[1] is x[a+10] ... and s[11] is x[a] for some a.  Carry contains
  // the last carry value (c[a+11]).
  //
  // The recursion relation (3) in Luscher's note says 
  //   delta[n] = x[n-s] = x[n-r] -c[n-1] or for n=a+12,
  //   delta[a+12] = x[a+7] - x[a] -c[a+11] where we use r=12, s=5 per eqn. (7)
  // This reduces to 
  // s[11] = s[4] - s[11] - carry.
  // The next number similarly will be given by s[10] = s[3] - s[10] - carry,
  // and so forth until s[0] is filled.
  // 
  // However, we need to skip 397, 202 or 109 numbers - these are not divisible 
  // by 12 - to "fare well in the spectral test".  

  advance(pDozens);

  // Since we wish at the end to have the 12 last numbers in the order of 
  // s[11] first, till s[0] last, we will have to do 1, 10, or 1 iterations 
  // and then re-arrange to place to get the oldest one in s[11].
  // Generically, this will imply re-arranging the s array at the end,
  // but we can treat the special case of endIters = 1 separately for superior
  // efficiency in the cases of levels 0 and 2.

  double  y1;

  if ( endIters == 1 ) {  	// Luxury levels 0 and 2 will go here
    y1 = randoms[ 4] - randoms[11] - carry;
    if ( y1 < 0.0 ) {
      y1 += 1.0;			
      carry = twoToMinus_48();
    } else {
      carry = 0.0;
    }
    randoms[11] = randoms[10];  
    randoms[10] = randoms[ 9];  
    randoms[ 9] = randoms[ 8];  
    randoms[ 8] = randoms[ 7];  
    randoms[ 7] = randoms[ 6];  
    randoms[ 6] = randoms[ 5];  
    randoms[ 5] = randoms[ 4];  
    randoms[ 4] = randoms[ 3];  
    randoms[ 3] = randoms[ 2];  
    randoms[ 2] = randoms[ 1];  
    randoms[ 1] = randoms[ 0];  
    randoms[ 0] = y1;

  } else {

    int m, nr, ns;
    for ( m = 0, nr = 11, ns = 4; m < endIters; ++m, --nr ) {
      y1 = randoms [ns] - randoms[nr] - carry;
      if ( y1 < 0.0 ) {
        y1 += 1.0;
        carry = twoToMinus_48();
      } else {
        carry = 0.0;
      }
      randoms[nr] = y1;
      --ns;
      if ( ns < 0 ) {
        ns = 11;
      }
    } // loop on m

    double temp[12];
    for (m=0; m<12; m++) {
      temp[m]=randoms[m];
    }

    ns = 11 - endIters;
    for (m=11; m>=0; --m) {
      randoms[m] = temp[ns];
      --ns;
      if ( ns < 0 ) {
        ns = 11;
      }
    } 

  }

  // Now when we return, there are 12 fresh usable numbers in s[11] ... s[0]

  index = 11;

} // update()

void Ranlux64Engine::advance(int dozens) {

  double  y1, y2, y3;
  double  cValue = twoToMinus_48();
  double  zero = 0.0;
  double  one  = 1.0;

		// Technical note:  We use Luscher's trick to only do the
		// carry subtraction when we really have to.  Like him, we use 
		// three registers instead of two so that we avoid sequences
		// like storing y1 then immediately replacing its value:
		// some architectures lose time when this is done.

  		// Luscher's ranlxd.c fills the stash going
		// upward.  We fill it downward to save a bit of time in the
		// flat() routine at no cost later.  This means that while
		// Luscher's ir is jr+5, our n-r is (n-s)-5.  (Note that
		// though ranlxd.c initializes ir and jr to 11 and 7, ir as
		// used is 5 more than jr because update is entered after 
		// incrementing ir.)  
		//

		// I have CAREFULLY checked that the algorithms do match
		// in all details.

  int k;
  for ( k = dozens; k > 0; --k ) {

    y1 = randoms[ 4] - randoms[11] - carry;

    y2 = randoms[ 3] - randoms[10];
    if ( y1 < zero ) {
      y1 += one;			
      y2 -= cValue;
    }
    randoms[11] = y1;

    y3 = randoms[ 2] - randoms[ 9];
    if ( y2 < zero ) {
      y2 += one;			
      y3 -= cValue;
    }
    randoms[10] = y2;

    y1 = randoms[ 1] - randoms[ 8];
    if ( y3 < zero ) {
      y3 += one;			
      y1 -= cValue;
    }
    randoms[ 9] = y3;

    y2 = randoms[ 0] - randoms[ 7];
    if ( y1 < zero ) {
      y1 += one;			
      y2 -= cValue;
    }
    randoms[ 8] = y1;

    y3 = randoms[11] - randoms[ 6];
    if ( y2 < zero ) {
      y2 += one;			
      y3 -= cValue;
    }
    randoms[ 7] = y2;

    y1 = randoms[10] - randoms[ 5];
    if ( y3 < zero ) {
      y3 += one;			
      y1 -= cValue;
    }
    randoms[ 6] = y3;

    y2 = randoms[ 9] - randoms[ 4];
    if ( y1 < zero ) {
      y1 += one;			
      y2 -= cValue;
    }
    randoms[ 5] = y1;

    y3 = randoms[ 8] - randoms[ 3];
    if ( y2 < zero ) {
      y2 += one;			
      y3 -= cValue;
    }
    randoms[ 4] = y2;

    y1 = randoms[ 7] - randoms[ 2];
    if ( y3 < zero ) {
      y3 += one;			
      y1 -= cValue;
    }
    randoms[ 3] = y3;

    y2 = randoms[ 6] - randoms[ 1];
    if ( y1 < zero ) {
      y1 += one;			
      y2 -= cValue;
    }
    randoms[ 2] = y1;

    y3 = randoms[ 5] - randoms[ 0];
    if ( y2 < zero ) {
      y2 += one;			
      y3 -= cValue;
    }
    randoms[ 1] = y2;

    if ( y3 < zero ) {
      y3 += one;			
      carry = cValue;
    }
    randoms[ 0] = y3; 

  } // End of major k loop doing 12 numbers at each cycle

} // advance(dozens)

void Ranlux64Engine::flatArray(const int size, double* vect) {
  for( int i=0; i < size; ++i ) {
    vect[i] = flat(); 
  }
}

void Ranlux64Engine::setSeed(long seed, int lux) {

// The initialization is carried out using a Multiplicative
// Congruential generator using formula constants of L'Ecuyer
// as described in "A review of pseudorandom number generators"
// (Fred James) published in Computer Physics Communications 60 (1990)
// pages 329-344

  const int ecuyer_a(53668);
  const int ecuyer_b(40014);
  const int ecuyer_c(12211);
  const int ecuyer_d(2147483563);

  const int lux_levels[3] = {109, 202, 397};
  theSeed = seed;

  if( (lux > 2)||(lux < 0) ){
     pDiscard = (lux >= 12) ? (lux-12) : lux_levels[1];
  }else{
     pDiscard = lux_levels[luxury];
  }
  pDozens  = pDiscard / 12;
  endIters = pDiscard % 12;

  long init_table[24];
  long next_seed = seed;
  long k_multiple;
  int i;
  next_seed &= 0xffffffff;
  while( next_seed >= ecuyer_d ) {
     next_seed -= ecuyer_d;
  }
  
  for(i = 0;i != 24;i++){
     k_multiple = next_seed / ecuyer_a;
     next_seed = ecuyer_b * (next_seed - k_multiple * ecuyer_a)
                                       - k_multiple * ecuyer_c;
     if(next_seed < 0) {
	next_seed += ecuyer_d;
     }
     next_seed &= 0xffffffff;
     init_table[i] = next_seed;
  } 
  // are we on a 64bit machine?
  if( sizeof(long) >= 8 ) {
     long topbits1, topbits2;
#ifdef WIN32
     topbits1 = ( seed >> 32) & 0xffff ;
     topbits2 = ( seed >> 48) & 0xffff ;
#else
     topbits1 = detail::rshift<32>(seed) & 0xffff ;
     topbits2 = detail::rshift<48>(seed) & 0xffff ;
#endif
     init_table[0] ^= topbits1;
     init_table[2] ^= topbits2;
     //std::cout << " init_table[0] " << init_table[0] << " from " << topbits1 << std::endl;
     //std::cout << " init_table[2] " << init_table[2] << " from " << topbits2 << std::endl;
  }   

  for(i = 0;i < 12; i++){
     randoms[i] = (init_table[2*i  ]      ) * 2.0 * twoToMinus_32() +
                  (init_table[2*i+1] >> 15) * twoToMinus_48();
     //if( randoms[i] < 0. || randoms[i]  > 1. ) {
     //std::cout << "setSeed:  init_table " << init_table[2*i  ] << std::endl;
     //std::cout << "setSeed:  init_table " << init_table[2*i+1] << std::endl;
     //std::cout << "setSeed:  random " << i << " is " << randoms[i] << std::endl;
     //}
  }

  carry = 0.0;
  if ( randoms[11] == 0. ) carry = twoToMinus_48();
  index = 11;

} // setSeed()

void Ranlux64Engine::setSeeds(const long * seeds, int lux) {
// old code only uses the first long in seeds
//  setSeed( *seeds ? *seeds : 32767, lux );
//  theSeeds = seeds;

// using code from Ranlux - even those are 32bit seeds, 
// that is good enough to completely differentiate the sequences

   const int ecuyer_a = 53668;
   const int ecuyer_b = 40014;
   const int ecuyer_c = 12211;
   const int ecuyer_d = 2147483563;

   const int lux_levels[3] = {109, 202, 397};
   const long *seedptr; 

   theSeeds = seeds;
   seedptr  = seeds;
 
   if(seeds == 0){
      setSeed(theSeed,lux);
      theSeeds = &theSeed;
      return;
   }

   theSeed = *seeds;

// number of additional random numbers that need to be 'thrown away'
// every 24 numbers is set using luxury level variable.

  if( (lux > 2)||(lux < 0) ){
     pDiscard = (lux >= 12) ? (lux-12) : lux_levels[1];
  }else{
     pDiscard = lux_levels[luxury];
  }
  pDozens  = pDiscard / 12;
  endIters = pDiscard % 12;

  long init_table[24];
  long next_seed = *seeds;
  long k_multiple;
  int i;
      
  for( i = 0;(i != 24)&&(*seedptr != 0);i++){
      init_table[i] =  *seedptr & 0xffffffff;
      seedptr++;
  }		       

  if(i != 24){
     next_seed = init_table[i-1];
     for(;i != 24;i++){
	k_multiple = next_seed / ecuyer_a;
	next_seed = ecuyer_b * (next_seed - k_multiple * ecuyer_a)
                                	  - k_multiple * ecuyer_c;
	if(next_seed < 0) {
	   next_seed += ecuyer_d;
	}
	next_seed &= 0xffffffff;
	init_table[i] = next_seed;
     }    
  }

  for(i = 0;i < 12; i++){
     randoms[i] = (init_table[2*i  ]      ) * 2.0 * twoToMinus_32() +
                  (init_table[2*i+1] >> 15) * twoToMinus_48();
  }

  carry = 0.0;
  if ( randoms[11] == 0. ) carry = twoToMinus_48();
  index = 11;

}

void Ranlux64Engine::saveStatus( const char filename[] ) const
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

void Ranlux64Engine::restoreStatus( const char filename[] )
{
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
     for (int i=0; i<12; ++i) {
       inFile >> randoms[i];
     }
     inFile >> carry; inFile >> index;
     inFile >> luxury; inFile >> pDiscard;
     pDozens  = pDiscard / 12;
     endIters = pDiscard % 12;
   }
}

void Ranlux64Engine::showStatus() const
{
   std::cout << std::endl;
   std::cout << "--------- Ranlux engine status ---------" << std::endl;
   std::cout << " Initial seed = " << theSeed << std::endl;
   std::cout << " randoms[] = ";
   for (int i=0; i<12; ++i) {
     std::cout << randoms[i] << std::endl;
   }
   std::cout << std::endl;
   std::cout << " carry = " << carry << ", index = " << index << std::endl;
   std::cout << " luxury = " << luxury << " pDiscard = " 
						<< pDiscard << std::endl;
   std::cout << "----------------------------------------" << std::endl;
}

std::ostream & Ranlux64Engine::put( std::ostream& os ) const
{
   char beginMarker[] = "Ranlux64Engine-begin";
  os << beginMarker << "\nUvec\n";
  std::vector<unsigned long> v = put();
  for (unsigned int i=0; i<v.size(); ++i) {
     os <<  v[i] <<  "\n";
  }
  return os;  
}

std::vector<unsigned long> Ranlux64Engine::put () const {
  std::vector<unsigned long> v;
  v.push_back (engineIDulong<Ranlux64Engine>());
  std::vector<unsigned long> t;
  for (int i=0; i<12; ++i) {
    t = DoubConv::dto2longs(randoms[i]);
    v.push_back(t[0]); v.push_back(t[1]);
  }
  t = DoubConv::dto2longs(carry);
  v.push_back(t[0]); v.push_back(t[1]);
  v.push_back(static_cast<unsigned long>(index));
  v.push_back(static_cast<unsigned long>(luxury));
  v.push_back(static_cast<unsigned long>(pDiscard));
  return v;
}

std::istream & Ranlux64Engine::get ( std::istream& is )
{
  char beginMarker [MarkerLen];
  is >> std::ws;
  is.width(MarkerLen);  // causes the next read to the char* to be <=
			// that many bytes, INCLUDING A TERMINATION \0 
			// (Stroustrup, section 21.3.2)
  is >> beginMarker;
  if (strcmp(beginMarker,"Ranlux64Engine-begin")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nInput stream mispositioned or"
	       << "\nRanlux64Engine state description missing or"
	       << "\nwrong engine type found." << std::endl;
     return is;
  }
  return getState(is);
}

std::string Ranlux64Engine::beginTag ( )  { 
  return "Ranlux64Engine-begin"; 
}

std::istream & Ranlux64Engine::getState ( std::istream& is )
{
  if ( possibleKeywordInput ( is, "Uvec", theSeed ) ) {
    std::vector<unsigned long> v;
    unsigned long uu;
    for (unsigned int ivec=0; ivec < VECTOR_STATE_SIZE; ++ivec) {
      is >> uu;
      if (!is) {
        is.clear(std::ios::badbit | is.rdstate());
        std::cerr << "\nRanlux64Engine state (vector) description improper."
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

  char endMarker   [MarkerLen];
  for (int i=0; i<12; ++i) {
     is >> randoms[i];
  }
  is >> carry; is >> index;
  is >> luxury; is >> pDiscard;
  pDozens  = pDiscard / 12;
  endIters = pDiscard % 12;
  is >> std::ws;
  is.width(MarkerLen);  
  is >> endMarker;
  if (strcmp(endMarker,"Ranlux64Engine-end")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nRanlux64Engine state description incomplete."
	       << "\nInput stream is probably mispositioned now." << std::endl;
     return is;
  }
  return is;
}

bool Ranlux64Engine::get (const std::vector<unsigned long> & v) {
  if ((v[0] & 0xffffffffUL) != engineIDulong<Ranlux64Engine>()) {
    std::cerr << 
    	"\nRanlux64Engine get:state vector has wrong ID word - state unchanged\n";
    return false;
  }
  return getState(v);
}

bool Ranlux64Engine::getState (const std::vector<unsigned long> & v) {
  if (v.size() != VECTOR_STATE_SIZE ) {
    std::cerr << 
    	"\nRanlux64Engine get:state vector has wrong length - state unchanged\n";
    return false;
  }
  std::vector<unsigned long> t(2);
  for (int i=0; i<12; ++i) {
    t[0] = v[2*i+1]; t[1] = v[2*i+2];
    randoms[i] = DoubConv::longs2double(t);
  }
  t[0] = v[25]; t[1] = v[26];
  carry    = DoubConv::longs2double(t);
  index    = v[27];
  luxury   = v[28];
  pDiscard = v[29]; 
  return true;
}

}  // namespace CLHEP
