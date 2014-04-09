// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                        --- RanluxEngine ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).
//
// Ranlux random number generator originally implemented in FORTRAN77
// by Fred James as part of the MATHLIB HEP library.
// 'RanluxEngine' is designed to fit into the CLHEP random number
// class structure.

// ===============================================================
// Adeyemi Adesanya - Created: 6th November 1995
// Gabriele Cosmo - Adapted & Revised: 22nd November 1995
// Adeyemi Adesanya - Added setSeeds() method: 2nd February 1996
// Gabriele Cosmo - Added flatArray() method: 8th February 1996
//                - Minor corrections: 31st October 1996
//                - Added methods for engine status: 19th November 1996
//                - Fixed bug in setSeeds(): 15th September 1997
// J.Marraffino   - Added stream operators and related constructor.
//                  Added automatic seed selection from seed table and
//                  engine counter: 14th Feb 1998
//                - Fixed bug: setSeeds() requires a zero terminated
//                  array of seeds: 19th Feb 1998
// Ken Smith      - Added conversion operators:  6th Aug 1998
// J. Marraffino  - Remove dependence on hepString class  13 May 1999
// M. Fischler    - In restore, checkFile for file not found    03 Dec 2004
// M. Fischler    - Methods put, getfor instance save/restore       12/8/04    
// M. Fischler    - split get() into tag validation and 
//                  getState() for anonymous restores           12/27/04    
// M. Fischler    - put/get for vectors of ulongs		3/14/05
// M. Fischler    - State-saving using only ints, for portability 4/12/05
//
// ===============================================================

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanluxEngine.h"
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

std::string RanluxEngine::name() const {return "RanluxEngine";}

RanluxEngine::RanluxEngine(long seed, int lux)
: HepRandomEngine()
{
   long seedlist[2]={0,0};

   luxury = lux;
   setSeed(seed, luxury);
   
   // setSeeds() wants a zero terminated array!
   seedlist[0]=theSeed;
   seedlist[1]=0;
   setSeeds(seedlist, luxury);
}

RanluxEngine::RanluxEngine()
: HepRandomEngine()
{
   long seed;
   long seedlist[2]={0,0};

   luxury = 3;
   int numEngines = numberOfEngines++;
   int cycle = std::abs(int(numEngines/maxIndex));
   int curIndex = std::abs(int(numEngines%maxIndex));

   long mask = ((cycle & 0x007fffff) << 8);
   HepRandom::getTheTableSeeds( seedlist, curIndex );
   seed = seedlist[0]^mask;
   setSeed(seed, luxury);
   
   // setSeeds() wants a zero terminated array!
   seedlist[0]=theSeed;
   seedlist[1]=0;
   setSeeds(seedlist, luxury);
}

RanluxEngine::RanluxEngine(int rowIndex, int colIndex, int lux)
: HepRandomEngine()
{
   long seed;
   long seedlist[2]={0,0};

   luxury = lux;
   int cycle = std::abs(int(rowIndex/maxIndex));
   int row = std::abs(int(rowIndex%maxIndex));
   int col = std::abs(int(colIndex%2));
   long mask = (( cycle & 0x000007ff ) << 20 );
   HepRandom::getTheTableSeeds( seedlist, row );
   seed = ( seedlist[col] )^mask;
   setSeed(seed, luxury);
   
   // setSeeds() wants a zero terminated array!
   seedlist[0]=theSeed;
   seedlist[1]=0;
   setSeeds(seedlist, luxury);
}

RanluxEngine::RanluxEngine( std::istream& is )
: HepRandomEngine()
{
  is >> *this;
}

RanluxEngine::~RanluxEngine() {}

void RanluxEngine::setSeed(long seed, int lux) {

// The initialisation is carried out using a Multiplicative
// Congruential generator using formula constants of L'Ecuyer 
// as described in "A review of pseudorandom number generators"
// (Fred James) published in Computer Physics Communications 60 (1990)
// pages 329-344

  const int ecuyer_a = 53668;
  const int ecuyer_b = 40014;
  const int ecuyer_c = 12211;
  const int ecuyer_d = 2147483563;

  const int lux_levels[5] = {0,24,73,199,365};  

  long int_seed_table[24];
  long next_seed = seed;
  long k_multiple;
  int i;
  
// number of additional random numbers that need to be 'thrown away'
// every 24 numbers is set using luxury level variable.

  theSeed = seed;
  if( (lux > 4)||(lux < 0) ){
     if(lux >= 24){
        nskip = lux - 24;
     }else{
        nskip = lux_levels[3]; // corresponds to default luxury level
     }
  }else{
     luxury = lux;
     nskip = lux_levels[luxury];
  }

   
  for(i = 0;i != 24;i++){
     k_multiple = next_seed / ecuyer_a;
     next_seed = ecuyer_b * (next_seed - k_multiple * ecuyer_a) 
     - k_multiple * ecuyer_c ;
     if(next_seed < 0)next_seed += ecuyer_d;
     int_seed_table[i] = next_seed % int_modulus;
  }     

  for(i = 0;i != 24;i++)
    float_seed_table[i] = int_seed_table[i] * mantissa_bit_24();

  i_lag = 23;
  j_lag = 9;
  carry = 0. ;

  if( float_seed_table[23] == 0. ) carry = mantissa_bit_24();
  
  count24 = 0;
}

void RanluxEngine::setSeeds(const long *seeds, int lux) {

   const int ecuyer_a = 53668;
   const int ecuyer_b = 40014;
   const int ecuyer_c = 12211;
   const int ecuyer_d = 2147483563;

   const int lux_levels[5] = {0,24,73,199,365};
   int i;
   long int_seed_table[24];
   long k_multiple,next_seed;
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

  if( (lux > 4)||(lux < 0) ){
     if(lux >= 24){
        nskip = lux - 24;
     }else{
        nskip = lux_levels[3]; // corresponds to default luxury level
     }
  }else{
     luxury = lux;
     nskip = lux_levels[luxury];
  }
      
  for( i = 0;(i != 24)&&(*seedptr != 0);i++){
      int_seed_table[i] = *seedptr % int_modulus;
      seedptr++;
  }		       

  if(i != 24){
     next_seed = int_seed_table[i-1];
     for(;i != 24;i++){
        k_multiple = next_seed / ecuyer_a;
        next_seed = ecuyer_b * (next_seed - k_multiple * ecuyer_a) 
        - k_multiple * ecuyer_c ;
        if(next_seed < 0)next_seed += ecuyer_d;
        int_seed_table[i] = next_seed % int_modulus;
     }          
  }

  for(i = 0;i != 24;i++)
    float_seed_table[i] = int_seed_table[i] * mantissa_bit_24();

  i_lag = 23;
  j_lag = 9;
  carry = 0. ;

  if( float_seed_table[23] == 0. ) carry = mantissa_bit_24();
  
  count24 = 0;
}

void RanluxEngine::saveStatus( const char filename[] ) const
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

void RanluxEngine::restoreStatus( const char filename[] )
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
        std::cerr << "\nRanluxEngine state (vector) description improper."
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
     for (int i=0; i<24; ++i)
       inFile >> float_seed_table[i];
     inFile >> i_lag; inFile >> j_lag;
     inFile >> carry; inFile >> count24;
     inFile >> luxury; inFile >> nskip;
   }
}

void RanluxEngine::showStatus() const
{
   std::cout << std::endl;
   std::cout << "--------- Ranlux engine status ---------" << std::endl;
   std::cout << " Initial seed = " << theSeed << std::endl;
   std::cout << " float_seed_table[] = ";
   for (int i=0; i<24; ++i)
     std::cout << float_seed_table[i] << " ";
   std::cout << std::endl;
   std::cout << " i_lag = " << i_lag << ", j_lag = " << j_lag << std::endl;
   std::cout << " carry = " << carry << ", count24 = " << count24 << std::endl;
   std::cout << " luxury = " << luxury << " nskip = " << nskip << std::endl;
   std::cout << "----------------------------------------" << std::endl;
}

double RanluxEngine::flat() {

  float next_random;
  float uni;
  int i;

  uni = float_seed_table[j_lag] - float_seed_table[i_lag] - carry;
  if(uni < 0. ){
     uni += 1.0;
     carry = mantissa_bit_24();
  }else{
     carry = 0.;
  }

  float_seed_table[i_lag] = uni;
  i_lag --;
  j_lag --;
  if(i_lag < 0) i_lag = 23;
  if(j_lag < 0) j_lag = 23;

  if( uni < mantissa_bit_12() ){
     uni += mantissa_bit_24() * float_seed_table[j_lag];
     if( uni == 0) uni = mantissa_bit_24() * mantissa_bit_24();
  }
  next_random = uni;
  count24 ++;

// every 24th number generation, several random numbers are generated
// and wasted depending upon the luxury level.

  if(count24 == 24 ){
     count24 = 0;
     for( i = 0; i != nskip ; i++){
         uni = float_seed_table[j_lag] - float_seed_table[i_lag] - carry;
         if(uni < 0. ){
            uni += 1.0;
            carry = mantissa_bit_24();
         }else{
            carry = 0.;
         }
         float_seed_table[i_lag] = uni;
	 i_lag --;
         j_lag --;
         if(i_lag < 0)i_lag = 23;
         if(j_lag < 0) j_lag = 23;
      }
  } 
  return (double) next_random;
}

void RanluxEngine::flatArray(const int size, double* vect)
{
  float next_random;
  float uni;
  int i;
  int index;

  for (index=0; index<size; ++index) {
    uni = float_seed_table[j_lag] - float_seed_table[i_lag] - carry;
    if(uni < 0. ){
       uni += 1.0;
       carry = mantissa_bit_24();
    }else{
       carry = 0.;
    }

    float_seed_table[i_lag] = uni;
    i_lag --;
    j_lag --;
    if(i_lag < 0) i_lag = 23;
    if(j_lag < 0) j_lag = 23;

    if( uni < mantissa_bit_12() ){
       uni += mantissa_bit_24() * float_seed_table[j_lag];
       if( uni == 0) uni = mantissa_bit_24() * mantissa_bit_24();
    }
    next_random = uni;
    vect[index] = (double)next_random;
    count24 ++;

// every 24th number generation, several random numbers are generated
// and wasted depending upon the luxury level.

    if(count24 == 24 ){
       count24 = 0;
       for( i = 0; i != nskip ; i++){
           uni = float_seed_table[j_lag] - float_seed_table[i_lag] - carry;
           if(uni < 0. ){
              uni += 1.0;
              carry = mantissa_bit_24();
           }else{
              carry = 0.;
           }
           float_seed_table[i_lag] = uni;
           i_lag --;
           j_lag --;
           if(i_lag < 0)i_lag = 23;
           if(j_lag < 0) j_lag = 23;
        }
    }
  }
} 

RanluxEngine::operator unsigned int() {
   return ((unsigned int)(flat() * exponent_bit_32()) & 0xffffffff) |
         (((unsigned int)(float_seed_table[i_lag]*exponent_bit_32())>>16) & 0xff);
   // needed because Ranlux doesn't fill all bits of the double
   // which therefore doesn't fill all bits of the integer.
}

std::ostream & RanluxEngine::put ( std::ostream& os ) const
{
   char beginMarker[] = "RanluxEngine-begin";
  os << beginMarker << "\nUvec\n";
  std::vector<unsigned long> v = put();
  for (unsigned int i=0; i<v.size(); ++i) {
     os <<  v[i] <<  "\n";
  }
  return os;  
}

std::vector<unsigned long> RanluxEngine::put () const {
  std::vector<unsigned long> v;
  v.push_back (engineIDulong<RanluxEngine>());
  for (int i=0; i<24; ++i) {
    v.push_back
    	(static_cast<unsigned long>(float_seed_table[i]/mantissa_bit_24()));
  }
  v.push_back(static_cast<unsigned long>(i_lag));
  v.push_back(static_cast<unsigned long>(j_lag));
  v.push_back(static_cast<unsigned long>(carry/mantissa_bit_24()));
  v.push_back(static_cast<unsigned long>(count24));
  v.push_back(static_cast<unsigned long>(luxury));
  v.push_back(static_cast<unsigned long>(nskip));
  return v;
}

std::istream & RanluxEngine::get ( std::istream& is )
{
  char beginMarker [MarkerLen];
  is >> std::ws;
  is.width(MarkerLen);  // causes the next read to the char* to be <=
			// that many bytes, INCLUDING A TERMINATION \0 
			// (Stroustrup, section 21.3.2)
  is >> beginMarker;
  if (strcmp(beginMarker,"RanluxEngine-begin")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nInput stream mispositioned or"
	       << "\nRanluxEngine state description missing or"
	       << "\nwrong engine type found." << std::endl;
     return is;
  }
  return getState(is);
}

std::string RanluxEngine::beginTag ( )  { 
  return "RanluxEngine-begin"; 
}

std::istream & RanluxEngine::getState ( std::istream& is )
{
  if ( possibleKeywordInput ( is, "Uvec", theSeed ) ) {
    std::vector<unsigned long> v;
    unsigned long uu;
    for (unsigned int ivec=0; ivec < VECTOR_STATE_SIZE; ++ivec) {
      is >> uu;
      if (!is) {
        is.clear(std::ios::badbit | is.rdstate());
        std::cerr << "\nRanluxEngine state (vector) description improper."
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
  for (int i=0; i<24; ++i) {
     is >> float_seed_table[i];
  }
  is >> i_lag; is >>  j_lag;
  is >> carry; is >> count24;
  is >> luxury; is >> nskip;
  is >> std::ws;
  is.width(MarkerLen);  
  is >> endMarker;
  if (strcmp(endMarker,"RanluxEngine-end")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nRanluxEngine state description incomplete."
	       << "\nInput stream is probably mispositioned now." << std::endl;
     return is;
  }
  return is;
}

bool RanluxEngine::get (const std::vector<unsigned long> & v) {
  if ((v[0] & 0xffffffffUL) != engineIDulong<RanluxEngine>()) {
    std::cerr << 
    	"\nRanluxEngine get:state vector has wrong ID word - state unchanged\n";
    return false;
  }
  return getState(v);
}

bool RanluxEngine::getState (const std::vector<unsigned long> & v) {
  if (v.size() != VECTOR_STATE_SIZE ) {
    std::cerr << 
    	"\nRanluxEngine get:state vector has wrong length - state unchanged\n";
    return false;
  }
  for (int i=0; i<24; ++i) {
    float_seed_table[i] = v[i+1]*mantissa_bit_24();
  }
  i_lag    = v[25];
  j_lag    = v[26];
  carry    = v[27]*mantissa_bit_24();
  count24  = v[28];
  luxury   = v[29];
  nskip    = v[30];
  return true;
}

}  // namespace CLHEP
