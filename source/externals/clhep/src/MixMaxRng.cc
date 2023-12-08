//
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                          HEP Random
//                       --- MixMaxRng ---
//                     class implementation file
// -----------------------------------------------------------------------
//
// This file interfaces the MixMax PseudoRandom Number Generator
// proposed by:
//
// G.K.Savvidy and N.G.Ter-Arutyunian,
//   On the Monte Carlo simulation of physical systems,
//   J.Comput.Phys. 97, 566 (1991);
//   Preprint EPI-865-16-86, Yerevan, Jan. 1986
//   http://dx.doi.org/10.1016/0021-9991(91)90015-D
//
// K.Savvidy
//   "The MIXMAX random number generator"
//   Comp. Phys. Commun. (2015)
//   http://dx.doi.org/10.1016/j.cpc.2015.06.003
//
// K.Savvidy and G.Savvidy
//   "Spectrum and Entropy of C-systems. MIXMAX random number generator"
//   Chaos, Solitons & Fractals, Volume 91, (2016) pp. 33-38
//   http://dx.doi.org/10.1016/j.chaos.2016.05.003
//
// =======================================================================
// Implementation by Konstantin Savvidy - Copyright 2004-2023
// July 2023 - Updated class structure upon suggestions from Marco Barbone
// September 2023 - fix (re-)initialization from Gabriele Cosmo
// =======================================================================

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/MixMaxRng.h"
#include "CLHEP/Random/engineIDulong.h"
#include "CLHEP/Utility/atomic_int.h"

#include <atomic>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string.h>        // for strcmp
#include <vector>

const unsigned long MASK32=0xffffffff;

namespace CLHEP {

namespace {
  // Number of instances with automatic seed selection
  CLHEP_ATOMIC_INT_TYPE numberOfEngines(0);
}

static const int MarkerLen = 64; // Enough room to hold a begin or end marker. 

MixMaxRng::MixMaxRng()
: HepRandomEngine()
{
   int numEngines = ++numberOfEngines;
   setSeed(static_cast<long>(numEngines));
}

MixMaxRng::MixMaxRng(long seed)
: HepRandomEngine()
{
   theSeed=seed;
   setSeed(seed);
}

MixMaxRng::MixMaxRng(std::istream& is)
: HepRandomEngine()
{
   get(is);
}

MixMaxRng::~MixMaxRng() 
{
}

MixMaxRng::MixMaxRng(const MixMaxRng& rng)
  : HepRandomEngine(rng)
{
   std::copy(rng.V, rng.V+N, V);
   sumtot= rng.sumtot;
   counter= rng.counter;
}

MixMaxRng& MixMaxRng::operator=(const MixMaxRng& rng)
{
   // Check assignment to self
   //
   if (this == &rng)  { return *this; }

   // Copy base class data
   //
   HepRandomEngine::operator=(rng);

   std::copy(rng.V, rng.V+N, V);
   sumtot= rng.sumtot;
   counter= rng.counter;

   return *this;
}

void MixMaxRng::saveStatus( const char filename[] ) const
{
   // Create a C file-handle or an object that can accept the same output
   FILE *fh= fopen(filename, "w");
   if( fh )
   {
     int j;
     fprintf(fh, "mixmax state, file version 1.0\n" );
     fprintf(fh, "N=%u; V[N]={", rng_get_N() );
     for (j=0; (j< (rng_get_N()-1) ); ++j) {
         fprintf(fh, "%llu, ", (unsigned long long)V[j] );
     }
     fprintf(fh, "%llu", (unsigned long long)V[rng_get_N()-1] );
     fprintf(fh, "}; " );
     fprintf(fh, "counter=%u; ", counter );
     fprintf(fh, "sumtot=%llu;\n", (unsigned long long)sumtot );
     fclose(fh);
   }
}

void MixMaxRng::restoreStatus( const char filename[] )
{
   // a function for reading the state from a file
   FILE* fin;
   if(  ( fin = fopen(filename, "r") ) )
   {
      char l=0;
      while ( l != '{' ) { // 0x7B = "{"
        l=fgetc(fin); // proceed until hitting opening bracket
      }
      ungetc(' ', fin);
   }
   else
   {
      fprintf(stderr, "mixmax -> read_state: error reading file %s\n", filename);
      throw std::runtime_error("Error in reading state file");
   }
    
   myuint_t vecVal;
   //printf("mixmax -> read_state: starting to read state from file\n");
   if (!fscanf(fin, "%llu", (unsigned long long*) &V[0]) )
   {
     fprintf(stderr, "mixmax -> read_state: error reading file %s\n", filename);
     throw std::runtime_error("Error in reading state file");
   }

   for (int i = 1; i < rng_get_N(); ++i)
   {
     if (!fscanf(fin, ", %llu", (unsigned long long*) &vecVal) )
     {
       fprintf(stderr, "mixmax -> read_state: error reading vector component i=%d from file %s\n", i, filename);
       throw std::runtime_error("Error in reading state file");
     }
     if( vecVal <= MixMaxRng::M61 )
     {
       V[i] = vecVal;
     }
     else
     {
       fprintf(stderr, "mixmax -> read_state: Invalid state vector value= %llu"
               " ( must be less than %llu ) "
               " obtained from reading file %s\n"
               , (unsigned long long)vecVal, (unsigned long long)MixMaxRng::M61, filename);
     }
   }
    
   int incounter;
   if (!fscanf( fin, "}; counter=%i; ", &incounter))
   {
     fprintf(stderr, "mixmax -> read_state: error reading counter from file %s\n", filename);
     throw std::runtime_error("Error in reading state file");
   }
   if( incounter <= rng_get_N() )
   {
     counter = incounter;
   }
   else
   {
     fprintf(stderr, "mixmax -> read_state: Invalid counter = %d"
             "  Must be 0 <= counter < %u\n" , counter, rng_get_N());
     print_state();
     throw std::runtime_error("Error in reading state counter");
   }
   precalc();
   myuint_t insumtot;
   if (!fscanf( fin, "sumtot=%llu\n", (unsigned long long*) &insumtot))
   {
     fprintf(stderr, "mixmax -> read_state: error reading checksum from file %s\n", filename);
     throw std::runtime_error("Error in reading state file");
   }

   if (sumtot != insumtot)
   {
     fprintf(stderr, "mixmax -> checksum error while reading state from file %s - corrupted?\n", filename);
     throw std::runtime_error("Error in reading state checksum");
   }
   fclose(fin);
}

void MixMaxRng::showStatus() const
{
   std::cout << std::endl;
   std::cout << "------- MixMaxRng engine status -------" << std::endl;

   std::cout << " Current state vector is:" << std::endl;
   print_state();
   std::cout << "---------------------------------------" << std::endl;
}

//  Preferred Seeding method
//  the values of 'Seeds' must be valid 32-bit integers 
//  Higher order bits will be ignored!!

void MixMaxRng::setSeeds(const long* Seeds, int seedNum)
{
   myID_t seed0, seed1= 0, seed2= 0, seed3= 0;

   if( seedNum < 1 ) {  // Assuming at least 2 seeds in vector...
       seed0= static_cast<myID_t>(Seeds[0]) & MASK32;
       seed1= static_cast<myID_t>(Seeds[1]) & MASK32;
   } 
   else
   {
     if( seedNum < 4 ) {
       seed0= static_cast<myID_t>(Seeds[0]) & MASK32;
       if( seedNum > 1){ seed1= static_cast<myID_t>(Seeds[1]) & MASK32; }
       if( seedNum > 2){ seed2= static_cast<myID_t>(Seeds[2]) & MASK32; }
     }
     if( seedNum >= 4 ) {
       seed0= static_cast<myID_t>(Seeds[0]) & MASK32;
       seed1= static_cast<myID_t>(Seeds[1]) & MASK32;
       seed2= static_cast<myID_t>(Seeds[2]) & MASK32;
       seed3= static_cast<myID_t>(Seeds[3]) & MASK32;
     }
   }
   theSeed = Seeds[0];
   theSeeds = Seeds;
   seed_uniquestream(seed3, seed2, seed1, seed0);
}

std::string MixMaxRng::engineName()
{
   return "MixMaxRng";
}

constexpr int MixMaxRng::rng_get_N()
{
   return N;
}

void MixMaxRng::flatArray(const int size, double* vect )
{
   for (int i=0; i<size; ++i) { vect[i] = flat(); }
}

std::ostream & MixMaxRng::put ( std::ostream& os ) const
{
   char beginMarker[] = "MixMaxRng-begin";
   char endMarker[]   = "MixMaxRng-end";

   long pr = os.precision(24);
   os << beginMarker << " ";
   os << theSeed << "\n";
   for (int i=0; i<rng_get_N(); ++i) {
      os <<  V[i] << "\n";
   }
   os << counter << "\n";
   os << sumtot << "\n";
   os << endMarker << "\n";
   os.precision(pr);
   return os;  
}

std::vector<unsigned long> MixMaxRng::put () const
{
   std::vector<unsigned long> vec;
   vec.push_back (engineIDulong<MixMaxRng>());
   for (int i=0; i<rng_get_N(); ++i)
   {
     vec.push_back(static_cast<unsigned long>(V[i] & MASK32));
       // little-ended order on all platforms
     vec.push_back(static_cast<unsigned long>(V[i] >> 32  ));
       // pack uint64 into a data structure which is 32-bit on some platforms
   }
   vec.push_back(static_cast<unsigned long>(counter));
   vec.push_back(static_cast<unsigned long>(sumtot & MASK32));
   vec.push_back(static_cast<unsigned long>(sumtot >> 32));
   return vec;
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
   for (int i=0; i<rng_get_N(); ++i)  is >> V[i];
   is >> counter;
   myuint_t checksum;
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
   if ( counter < 0 || counter > rng_get_N() ) {
       std::cerr << "\nMixMaxRng::getState(): "
                 << "vector read wrong value of counter from file!"
                 << "\nInput stream is probably mispositioned now.\n";
       return is;
   }
   precalc();
   if ( checksum != sumtot) {
       std::cerr << "\nMixMaxRng::getState(): "
                 << "checksum disagrees with value stored in file!"
                 << "\nInput stream is probably mispositioned now.\n";
       return is;
   }
   return is;
}

bool MixMaxRng::get (const std::vector<unsigned long> & vec)
{
   if ((vec[0] & 0xffffffffUL) != engineIDulong<MixMaxRng>())
   {
     std::cerr << 
        "\nMixMaxRng::get(): vector has wrong ID word - state unchanged\n";
     return false;
   }
   return getState(vec);
}

bool MixMaxRng::getState (const std::vector<unsigned long> & vec)
{
   if (vec.size() != VECTOR_STATE_SIZE ) {
     std::cerr <<
        "\nMixMaxRng::getState(): vector has wrong length - state unchanged\n";
     return false;
   }
   for (int i=1; i<2*rng_get_N() ; i=i+2) {
     V[i/2]= ( (vec[i] & MASK32) | ( (myuint_t)(vec[i+1]) << 32 ) );
     // unpack from a data structure which is 32-bit on some platforms
   }
   counter = (int)vec[2*rng_get_N()+1];
   precalc();
   if ( ( (vec[2*rng_get_N()+2] & MASK32)
        | ( (myuint_t)(vec[2*rng_get_N()+3]) << 32 ) ) != sumtot) {
     std::cerr << "\nMixMaxRng::getState(): vector has wrong checksum!"
               << "\nInput vector is probably mispositioned now.\n";
     return false;
   }
   return true;
}

myuint_t MixMaxRng::iterate_raw_vec(myuint_t* Y, myuint_t sumtotOld)
{
   // operates with a raw vector, uses known sum of elements of Y
            
   myuint_t  tempP, tempV;
   Y[0] = ( tempV = sumtotOld);
   myuint_t insumtot = Y[0], ovflow = 0; // will keep a running sum of all new elements
   tempP = 0;              // will keep a partial sum of all old elements
   for (int i=1; (i<N); ++i)
   {
     myuint_t tempPO = MULWU(tempP);
     tempP = modadd(tempP, Y[i]);
     tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); // new Y[i] = old Y[i] + old partial * m
     Y[i] = tempV;
     insumtot += tempV; if (insumtot < tempV) {++ovflow;}
   }
   return MIXMAX_MOD_MERSENNE(MIXMAX_MOD_MERSENNE(insumtot) + (ovflow <<3 ));
}
        
myuint_t MixMaxRng::get_next()
{
   int i = counter;
            
   if ((i<=(N-1)) )
   {
     ++counter;
     return V[i];
   }
   else
   {
     sumtot = iterate_raw_vec(V, sumtot);
     counter=2;
     return V[1];
   }
}

myuint_t MixMaxRng::precalc()
{
   myuint_t temp = 0;
   for (int i=0; i < N; ++i) { temp = MIXMAX_MOD_MERSENNE(temp + V[i]); }
   sumtot = temp;
   return temp;
}
        
void MixMaxRng::state_init()
{
   for (int i=1; i < N; ++i) { V[i] = 0; }
   V[0] = 1;
   counter = N;  // set the counter to N if iteration should happen right away
   sumtot = 1;
}

void MixMaxRng::seed_spbox(myuint_t seed)
{
   // a 64-bit LCG from Knuth line 26, in combination with a bit swap is used to seed

   if (seed == 0)
     throw std::runtime_error("try seeding with nonzero seed next time");

   const myuint_t MULT64=6364136223846793005ULL;
   sumtot = 0;
            
   myuint_t l = seed;
            
   for (int i=0; i < N; ++i)
   {
     l*=MULT64; l = (l << 32) ^ (l>>32);
     V[i] = l & M61;
     sumtot = MIXMAX_MOD_MERSENNE(sumtot + V[(i)]); 
   }
   counter = N;  // set the counter to N if iteration should happen right away
}

void MixMaxRng::seed_uniquestream( myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID )
{
   state_init();
   sumtot = apply_bigskip(V, V, clusterID, machineID, runID, streamID );
   counter = 1;
}
        
myuint_t MixMaxRng::apply_bigskip( myuint_t* Vout, myuint_t* Vin, myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID )
{
  /*
    makes a derived state vector, Vout, from the mother state vector Vin
    by skipping a large number of steps, determined by the given seeding ID's
             
    it is mathematically guaranteed that the substreams derived in this way from the SAME (!!!) Vin will not collide provided
    1) at least one bit of ID is different
    2) less than 10^100 numbers are drawn from the stream
    (this is good enough : a single CPU will not exceed this in the lifetime of the universe, 10^19 sec,
    even if it had a clock cycle of Planch time, 10^44 Hz )
             
    Caution: never apply this to a derived vector, just choose some mother vector Vin, for example the unit vector by state_init(),
    and use it in all your runs, just change runID to get completely nonoverlapping streams of random numbers on a different day.
             
    clusterID and machineID are provided for the benefit of large organizations who wish to ensure that a simulation
    which is running in parallel on a large number of  clusters and machines will have non-colliding source of random numbers.
             
    did i repeat it enough times? the non-collision guarantee is absolute, not probabilistic

   */

   const myuint_t skipMat17[128][17] =
   #include "CLHEP/Random/mixmax_skip_N17.icc"
   ;
            
   const myuint_t* skipMat[128];
   for (int i=0; i<128; ++i) { skipMat[i] = skipMat17[i]; }
            
   myID_t IDvec[4] = {streamID, runID, machineID, clusterID};
   int r,i,j,  IDindex;
   myID_t id;
   myuint_t Y[N], cum[N];
   myuint_t coeff;
   myuint_t* rowPtr;
   myuint_t insumtot=0;
            
   for (i=0; i<N; ++i) { Y[i] = Vin[i]; insumtot = modadd( insumtot, Vin[i]); }
   for (IDindex=0; IDindex<4; ++IDindex)
   { // go from lower order to higher order ID
     id=IDvec[IDindex];
     //printf("now doing ID at level %d, with ID = %d\n", IDindex, id);
     r = 0;
     while (id)
     {
       if (id & 1)
       {
         rowPtr = (myuint_t*)skipMat[r + IDindex*8*sizeof(myID_t)];
         for (i=0; i<N; ++i){ cum[i] = 0; }
         for (j=0; j<N; ++j)
         { // j is lag, enumerates terms of the poly
           // for zero lag Y is already given
           coeff = rowPtr[j]; // same coeff for all i
           for (i =0; i<N; ++i){
             cum[i] =  fmodmulM61( cum[i], coeff ,  Y[i] ) ;
           }
           insumtot = iterate_raw_vec(Y, insumtot);
         }
         insumtot=0;
         for (i=0; i<N; ++i){ Y[i] = cum[i]; insumtot = modadd( insumtot, cum[i]); } ;
       }
       id = (id >> 1); ++r; // bring up the r-th bit in the ID
     }
   }
   insumtot=0;
   for (i=0; i<N; ++i){ Vout[i] = Y[i]; insumtot = modadd( insumtot, Y[i]); }
   // returns sumtot, and copy the vector over to Vout
   return insumtot;
}
        
#if defined(__x86_64__)
  myuint_t MixMaxRng::mod128(__uint128_t s)
  {
    myuint_t s1;
    s1 = ( (  ((myuint_t)s)&M61 ) + (  ((myuint_t)(s>>64)) * 8 ) + ( ((myuint_t)s) >>BITS) );
    return	MIXMAX_MOD_MERSENNE(s1);
  }
  myuint_t MixMaxRng::fmodmulM61(myuint_t cum, myuint_t a, myuint_t b)
  {
    __uint128_t temp;
    temp = (__uint128_t)a*(__uint128_t)b + cum;
    return mod128(temp);
  }
#else // on all other platforms, including 32-bit linux, PPC and PPC64, ARM and all Windows
  myuint_t MixMaxRng::fmodmulM61(myuint_t cum, myuint_t s, myuint_t a)
  {
    const myuint_t MASK32=0xFFFFFFFFULL;
    myuint_t o,ph,pl,ah,al;
    o=(s)*a;
    ph = ((s)>>32);
    pl = (s) & MASK32;
    ah = a>>32;
    al = a & MASK32;
    o = (o & M61) + ((ph*ah)<<3) + ((ah*pl+al*ph + ((al*pl)>>32))>>29) ;
    o += cum;
    o = (o & M61) + ((o>>61));
    return o;
  }
#endif

void MixMaxRng::print_state() const
{
   std::cout << "mixmax state, file version 1.0\n";
   std::cout << "N=" << rng_get_N() << "; V[N]={";
   for (int j=0; (j< (rng_get_N()-1) ); ++j) { std::cout << V[j] << ", "; }
   std::cout << V[rng_get_N()-1];
   std::cout << "}; ";
   std::cout << "counter= " << counter;
   std::cout << "sumtot= " << sumtot << "\n";
}

MixMaxRng MixMaxRng::Branch()
{
   sumtot = iterate_raw_vec(V, sumtot); counter = 1;
   MixMaxRng tmp=*this;
   tmp.BranchInplace(0); // daughter id
   return tmp;
}
    
void MixMaxRng::BranchInplace(int id)
{
   // Dont forget to iterate the mother, when branching the daughter, or else will have collisions!
   // a 64-bit LCG from Knuth line 26, is used to mangle a vector component
   constexpr myuint_t MULT64=6364136223846793005ULL;
   myuint_t tmp=V[id];
   V[1] *= MULT64; V[id] &= M61;
   sumtot = MIXMAX_MOD_MERSENNE( sumtot + V[id] - tmp + M61);
   sumtot = iterate_raw_vec(V, sumtot);
   counter = 1;
}

}  // namespace CLHEP
