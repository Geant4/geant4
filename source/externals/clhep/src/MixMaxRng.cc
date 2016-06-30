// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                          HEP Random
//                       --- MixMaxRng ---
//                     class implementation file
// -----------------------------------------------------------------------
//
// This file interfaces the PseudoRandom Number Generator 
// proposed by  N.Z. Akopov, G.K.Saviddy & N.G. Ter-Arutyunian
//   "Matrix Generator of Pseudorandom Numbers"
//       J. Compt. Phy. 97, 573 (1991)
//  Preprint: EPI-867(18)-86, Yerevan June 1986.
//
// Implementation by Konstantin Savvidy - 2004-2015
//        "The MIXMAX random number generator"
//        Comp. Phys. Commun. (2015)
//        http://dx.doi.org/10.1016/j.cpc.2015.06.003
//
//  Release 0.99 and later: released under the LGPL license version 3.0
// =======================================================================
// CLHEP interface implemented by 
//     J. Apostolakis, G. Cosmo & K. Savvidy - Created: 6th July 2015
//     CLHEP interface released under the LGPL license version 3.0
// =======================================================================

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/MixMaxRng.h"
#include "CLHEP/Random/engineIDulong.h"
#include "CLHEP/Utility/atomic_int.h"

#include <string.h>        // for strcmp
#include <cmath>
#include <cstdlib>
#include <stdint.h>

#include "CLHEP/Random/mixmax.h"

const unsigned long MASK32=0xffffffff;

namespace CLHEP {

namespace {
  // Number of instances with automatic seed selection
  CLHEP_ATOMIC_INT_TYPE numberOfEngines(0);
}

static const int MarkerLen = 64; // Enough room to hold a begin or end marker. 

std::string MixMaxRng::name() const { return "MixMaxRng"; } // N=" + N

MixMaxRng::MixMaxRng()
: HepRandomEngine()
{
   int numEngines = ++numberOfEngines;
   fRngState= rng_alloc();
   setSeed(static_cast<long>(numEngines));
}

MixMaxRng::MixMaxRng(long seed)
: HepRandomEngine()
{
   fRngState= rng_alloc();
   setSeed(seed);
}

MixMaxRng::MixMaxRng(std::istream& is)
: HepRandomEngine()
{
   get(is);
}

MixMaxRng::~MixMaxRng() 
{
   rng_free( fRngState ); 
}

MixMaxRng::MixMaxRng(const MixMaxRng& rng)
  : HepRandomEngine(rng)
{
   fRngState= rng_copy( rng.fRngState->V );
   fRngState->sumtot= rng.fRngState->sumtot;
   fRngState->counter= rng.fRngState->counter;
}

MixMaxRng& MixMaxRng::operator=(const MixMaxRng& rng)
{
   // Check assignment to self
   //
   if (this == &rng)  { return *this; }

   // Copy base class data
   //
   HepRandomEngine::operator=(rng);

   // Copy data
   //
   rng_free( fRngState );
   fRngState= rng_copy( rng.fRngState->V );
   fRngState->sumtot= rng.fRngState->sumtot;
   fRngState->counter= rng.fRngState->counter;

   return *this;
}

void MixMaxRng::saveStatus( const char filename[] ) const
{
   // Create a C file-handle or an object that can accept the same output
   FILE *fh= fopen(filename, "w");
   if( fh )
   {
     fRngState->fh= fh;
     print_state(fRngState);
     fclose(fh);
   }
   fRngState->fh= 0;
}

void MixMaxRng::restoreStatus( const char filename[] )
{
   read_state(fRngState, filename);
}

void MixMaxRng::showStatus() const
{
   std::cout << std::endl;
   std::cout << "------- MixMaxRng engine status -------" << std::endl;

   std::cout << " Current state vector is:" << std::endl;
   fRngState->fh=stdout;
   print_state(fRngState);
   std::cout << "---------------------------------------" << std::endl;
}

void MixMaxRng::setSeed(long longSeed, int /* extraSeed */)
{
   unsigned long seed0;

   theSeed = longSeed;
   if( sizeof(long) > 4)  // C standard says long is at least 32-bits
     seed0= static_cast<unsigned long>(longSeed) & MASK32 ;
   else
     seed0= longSeed;

   seed_spbox(fRngState, seed0);
}

//  Preferred Seeding method
//  the values of 'Seeds' must be valid 32-bit integers 
//  Higher order bits will be ignored!!

void MixMaxRng::setSeeds(const long* Seeds, int seedNum)
{
   unsigned long seed0, seed1= 0, seed2= 0, seed3= 0;

   if( seedNum < 1 ) {  // Assuming at least 2 seeds in vector...
       seed0= static_cast<unsigned long>(Seeds[0]) & MASK32;
       seed1= static_cast<unsigned long>(Seeds[1]) & MASK32;
   } 
   else
   {
     if( seedNum < 4 ) {
       seed0= static_cast<unsigned long>(Seeds[0]) & MASK32;
       if( seedNum > 1){ seed1= static_cast<unsigned long>(Seeds[1]) & MASK32; }
       if( seedNum > 2){ seed2= static_cast<unsigned long>(Seeds[2]) & MASK32; }
     }
     if( seedNum >= 4 ) {
       seed0= static_cast<unsigned long>(Seeds[0]) & MASK32;
       seed1= static_cast<unsigned long>(Seeds[1]) & MASK32;
       seed2= static_cast<unsigned long>(Seeds[2]) & MASK32;
       seed3= static_cast<unsigned long>(Seeds[3]) & MASK32;
     }
   }
   theSeed = Seeds[0];
   theSeeds = Seeds;
   seed_uniquestream(fRngState, seed3, seed2, seed1, seed0);
}

double MixMaxRng::flat()
{
   return get_next_float(fRngState);
}

void MixMaxRng::flatArray(const int size, double* vect )
{
   // fill_array( fRngState, size, arrayDbl );
   for (int i=0; i<size; ++i) { vect[i] = flat(); }
}

MixMaxRng::operator unsigned int()
{
   return static_cast<unsigned int>(get_next(fRngState));
   // get_next returns a 64-bit integer, of which the lower 61 bits
   // are random and upper 3 bits are zero
}

std::ostream & MixMaxRng::put ( std::ostream& os ) const
{
   char beginMarker[] = "MixMaxRng-begin";
   char endMarker[]   = "MixMaxRng-end";

   int pr = os.precision(24);
   os << beginMarker << " ";
   os << theSeed << " ";
   for (int i=0; i<rng_get_N(); ++i) {
      os <<  fRngState->V[i] << "\n";
   }
   os << fRngState->counter << "\n";
   os << fRngState->sumtot << "\n";
   os << endMarker << "\n";
   os.precision(pr);
   return os;  
}

std::vector<unsigned long> MixMaxRng::put () const
{
   std::vector<unsigned long> v;
   v.push_back (engineIDulong<MixMaxRng>());
   for (int i=0; i<rng_get_N(); ++i) {
     v.push_back(static_cast<unsigned long>(fRngState->V[i] & MASK32));
       // little-ended order on all platforms
     v.push_back(static_cast<unsigned long>(fRngState->V[i] >> 32  ));
       // pack uint64 into a data structure which is 32-bit on some platforms
   }
   v.push_back(static_cast<unsigned long>(fRngState->counter));
   v.push_back(static_cast<unsigned long>(fRngState->sumtot & MASK32));
   v.push_back(static_cast<unsigned long>(fRngState->sumtot >> 32));
   return v;
}

std::istream & MixMaxRng::get  ( std::istream& is)
{
   char beginMarker [MarkerLen];
   is >> std::ws;
   is.width(MarkerLen);  // causes the next read to the char* to be <=
                         // that many bytes, INCLUDING A TERMINATION \0 
                         // (Stroustrup, section 21.3.2)
   is >> beginMarker;
   if (strcmp(beginMarker,"MixMaxRng-begin")) {
      is.clear(std::ios::badbit | is.rdstate());
      std::cerr << "\nInput stream mispositioned or"
                << "\nMixMaxRng state description missing or"
                << "\nwrong engine type found." << std::endl;
      return is;
   }
   return getState(is);
}

std::string MixMaxRng::beginTag ()
{ 
   return "MixMaxRng-begin"; 
}

std::istream &  MixMaxRng::getState ( std::istream& is )
{
   char endMarker[MarkerLen];
   is >> theSeed;
   for (int i=0; i<rng_get_N(); ++i)  is >> fRngState->V[i];
   is >> fRngState->counter;
   myuint checksum;
   is >> checksum;
   is >> std::ws;
   is.width(MarkerLen);
   is >> endMarker;
   if (strcmp(endMarker,"MixMaxRng-end")) {
       is.clear(std::ios::badbit | is.rdstate());
       std::cerr << "\nMixMaxRng state description incomplete."
                 << "\nInput stream is probably mispositioned now.\n";
       return is;
   }
   if ( fRngState->counter < 0 || fRngState->counter > rng_get_N() ) {
       std::cerr << "\nMixMaxRng::getState(): "
                 << "vector read wrong value of counter from file!"
                 << "\nInput stream is probably mispositioned now.\n";
       return is;
   }
   precalc(fRngState);
   if ( checksum != fRngState->sumtot) {
       std::cerr << "\nMixMaxRng::getState(): "
                 << "checksum disagrees with value stored in file!"
                 << "\nInput stream is probably mispositioned now.\n";
       return is;
   }
   return is;
}

bool MixMaxRng::get (const std::vector<unsigned long> & v)
{
   if ((v[0] & 0xffffffffUL) != engineIDulong<MixMaxRng>()) {
     std::cerr << 
        "\nMixMaxRng::get(): vector has wrong ID word - state unchanged\n";
     return false;
   }
   return getState(v);
}

bool MixMaxRng::getState (const std::vector<unsigned long> & v)
{
   if (v.size() != VECTOR_STATE_SIZE ) {
     std::cerr <<
        "\nMixMaxRng::getState(): vector has wrong length - state unchanged\n";
     return false;
   }
   for (int i=1; i<2*rng_get_N() ; i=i+2) {
     fRngState->V[i/2]= ( (v[i] & MASK32) | ( (myuint)(v[i+1]) << 32 ) );
     // unpack from a data structure which is 32-bit on some platforms
   }
   fRngState->counter = v[2*rng_get_N()+1];
   precalc(fRngState);
   if ( ( (v[2*rng_get_N()+2] & MASK32)
        | ( (myuint)(v[2*rng_get_N()+3]) << 32 ) ) != fRngState->sumtot) {
     std::cerr << "\nMixMaxRng::getState(): vector has wrong checksum!"
               << "\nInput vector is probably mispositioned now.\n";
     return false;
   }
   return true;
}

}  // namespace CLHEP
