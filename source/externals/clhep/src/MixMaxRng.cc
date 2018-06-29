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
// Implementation by Konstantin Savvidy - Copyright 2004-2017
// =======================================================================

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/MixMaxRng.h"
#include "CLHEP/Random/engineIDulong.h"
#include "CLHEP/Utility/atomic_int.h"

#include <string.h>        // for strcmp
#include <cmath>

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
   S.V = rng.S.V;
   S.sumtot= rng.S.sumtot;
   S.counter= rng.S.counter;
}

MixMaxRng& MixMaxRng::operator=(const MixMaxRng& rng)
{
   // Check assignment to self
   //
   if (this == &rng)  { return *this; }

   // Copy base class data
   //
   HepRandomEngine::operator=(rng);

   S.V = rng.S.V;
   S.sumtot= rng.S.sumtot;
   S.counter= rng.S.counter;

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
     for (j=0; (j< (rng_get_N()-1) ); j++) {
         fprintf(fh, "%llu, ", S.V[j] );
     }
     fprintf(fh, "%llu", S.V[rng_get_N()-1] );
     fprintf(fh, "}; " );
     fprintf(fh, "counter=%u; ", S.counter );
     fprintf(fh, "sumtot=%llu;\n", S.sumtot );
     fclose(fh);
   }
}

#define MIXMAX_ARRAY_INDEX_OUT_OF_BOUNDS   0xFF01
#define MIXMAX_SEED_WAS_ZERO               0xFF02
#define MIXMAX_ERROR_READING_STATE_FILE    0xFF03
#define MIXMAX_ERROR_READING_STATE_COUNTER       0xFF04
#define MIXMAX_ERROR_READING_STATE_CHECKSUM      0xFF05

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
      exit(MIXMAX_ERROR_READING_STATE_FILE);
   }
    
   myuint_t vecVal;
   //printf("mixmax -> read_state: starting to read state from file\n");
   if (!fscanf(fin, "%llu", &S.V[0]) )
   {
     fprintf(stderr, "mixmax -> read_state: error reading file %s\n", filename);
     exit(MIXMAX_ERROR_READING_STATE_FILE);
   }

   int i;
   for( i = 1; i < rng_get_N(); i++)
   {
     if (!fscanf(fin, ", %llu", &vecVal) )
     {
       fprintf(stderr, "mixmax -> read_state: error reading vector component i=%d from file %s\n", i, filename);
       exit(MIXMAX_ERROR_READING_STATE_FILE);
     }
     if(  vecVal <= MixMaxRng::M61 )
     {
       S.V[i] = vecVal;
     }
     else
     {
       fprintf(stderr, "mixmax -> read_state: Invalid state vector value= %llu"
               " ( must be less than %llu ) "
               " obtained from reading file %s\n"
               , vecVal, MixMaxRng::M61, filename);
     }
   }
    
   int counter;
   if (!fscanf( fin, "}; counter=%i; ", &counter))
   {
     fprintf(stderr, "mixmax -> read_state: error reading counter from file %s\n", filename);
     exit(MIXMAX_ERROR_READING_STATE_FILE);
   }
   if( counter <= rng_get_N() )
   {
     S.counter= counter;
   }
   else
   {
     fprintf(stderr, "mixmax -> read_state: Invalid counter = %d"
             "  Must be 0 <= counter < %u\n" , counter, rng_get_N());
     print_state();
     exit(MIXMAX_ERROR_READING_STATE_COUNTER);
   }
   precalc();
   myuint_t sumtot;
   if (!fscanf( fin, "sumtot=%llu\n", &sumtot))
   {
     fprintf(stderr, "mixmax -> read_state: error reading checksum from file %s\n", filename);
     exit(MIXMAX_ERROR_READING_STATE_FILE);
   }

   if (S.sumtot != sumtot)
   {
     fprintf(stderr, "mixmax -> checksum error while reading state from file %s - corrupted?\n", filename);
     exit(MIXMAX_ERROR_READING_STATE_CHECKSUM);
   }
   fclose(fin);
}

#undef MIXMAX_ARRAY_INDEX_OUT_OF_BOUNDS
#undef MIXMAX_SEED_WAS_ZERO
#undef MIXMAX_ERROR_READING_STATE_FILE
#undef MIXMAX_ERROR_READING_STATE_COUNTER
#undef MIXMAX_ERROR_READING_STATE_CHECKSUM

void MixMaxRng::showStatus() const
{
   std::cout << std::endl;
   std::cout << "------- MixMaxRng engine status -------" << std::endl;

   std::cout << " Current state vector is:" << std::endl;
   print_state();
   std::cout << "---------------------------------------" << std::endl;
}

void MixMaxRng::setSeed(long longSeed, int /* extraSeed */)
{
   //seed_uniquestream(0,0,0,longSeed);
   theSeed = longSeed;
   seed_spbox(longSeed);
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

constexpr long long int MixMaxRng::rng_get_SPECIAL()
{
   return SPECIAL;
}

constexpr int MixMaxRng::rng_get_SPECIALMUL()
{
   return SPECIALMUL;
}

double MixMaxRng::generate(int i)
{
   S.counter++;
#if defined(__clang__) || defined(__llvm__)
   return INV_M61*static_cast<double>(S.V[i]);
#elif defined(__GNUC__) && (__GNUC__ < 7) && (!defined(__ICC)) && defined(__x86_64__) && defined(__SSE2_MATH__)
   int64_t Z=S.V[i];
   double F=0.0;
   //#warning Using the inline assembler
    /* using SSE inline assemly to zero the xmm register, just before int64 -> double conversion,
       not necessary in GCC-5 or better, but huge penalty on earlier compilers
    */
   __asm__ __volatile__(  "pxor %0, %0;"
                          "cvtsi2sdq %1, %0;"
                          :"=x"(F)
                          :"r"(Z)
                       );
   return F*INV_M61;
#else
  //#warning other method
   return convert1double(S.V[i]); //get_next_float_packbits();
#endif
}

double MixMaxRng::iterate()
{
   myuint_t* Y=S.V.data();
   myuint_t  tempP, tempV;
   Y[0] = ( tempV = S.sumtot);
   myuint_t sumtot = Y[0], ovflow = 0; // will keep a running sum of all new elements
   tempP = 0;              // will keep a partial sum of all old elements
   myuint_t tempPO;
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[1] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[1]  = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[2] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[2]  = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[3] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[3]  = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[4] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[4]  = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[5] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[5]  = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[6] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[6]  = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[7] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[7]  = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[8] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[8]  = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[9] ); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[9]  = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[10]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[10] = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[11]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[11] = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[12]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[12] = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[13]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[13] = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[14]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[14] = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[15]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[15] = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   tempPO = MULWU(tempP); tempP = modadd(tempP, Y[16]); tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); Y[16] = tempV; sumtot += tempV; if (sumtot < tempV) {ovflow++;};
   S.sumtot = MIXMAX_MOD_MERSENNE(MIXMAX_MOD_MERSENNE(sumtot) + (ovflow <<3 ));

   S.counter=2;
   return double(S.V[1])*INV_M61;
}

void MixMaxRng::flatArray(const int size, double* vect )
{
   // fill_array( S, size, arrayDbl );
   for (int i=0; i<size; ++i) { vect[i] = flat(); }
}

MixMaxRng::operator double()
{
  return flat();
}

MixMaxRng::operator float()
{
  return float( flat() );
}

MixMaxRng::operator unsigned int()
{
   return static_cast<unsigned int>(get_next());
   // clhep_get_next returns a 64-bit integer, of which the lower 61 bits
   // are random and upper 3 bits are zero
}

std::ostream & MixMaxRng::put ( std::ostream& os ) const
{
   char beginMarker[] = "MixMaxRng-begin";
   char endMarker[]   = "MixMaxRng-end";

   int pr = os.precision(24);
   os << beginMarker << " ";
   os << theSeed << "\n";
   for (int i=0; i<rng_get_N(); ++i) {
      os <<  S.V[i] << "\n";
   }
   os << S.counter << "\n";
   os << S.sumtot << "\n";
   os << endMarker << "\n";
   os.precision(pr);
   return os;  
}

std::vector<unsigned long> MixMaxRng::put () const
{
   std::vector<unsigned long> v;
   v.push_back (engineIDulong<MixMaxRng>());
   for (int i=0; i<rng_get_N(); ++i)
   {
     v.push_back(static_cast<unsigned long>(S.V[i] & MASK32));
       // little-ended order on all platforms
     v.push_back(static_cast<unsigned long>(S.V[i] >> 32  ));
       // pack uint64 into a data structure which is 32-bit on some platforms
   }
   v.push_back(static_cast<unsigned long>(S.counter));
   v.push_back(static_cast<unsigned long>(S.sumtot & MASK32));
   v.push_back(static_cast<unsigned long>(S.sumtot >> 32));
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
   for (int i=0; i<rng_get_N(); ++i)  is >> S.V[i];
   is >> S.counter;
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
   if ( S.counter < 0 || S.counter > rng_get_N() ) {
       std::cerr << "\nMixMaxRng::getState(): "
                 << "vector read wrong value of counter from file!"
                 << "\nInput stream is probably mispositioned now.\n";
       return is;
   }
   precalc();
   if ( checksum != S.sumtot) {
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
     S.V[i/2]= ( (v[i] & MASK32) | ( (myuint_t)(v[i+1]) << 32 ) );
     // unpack from a data structure which is 32-bit on some platforms
   }
   S.counter = v[2*rng_get_N()+1];
   precalc();
   if ( ( (v[2*rng_get_N()+2] & MASK32)
        | ( (myuint_t)(v[2*rng_get_N()+3]) << 32 ) ) != S.sumtot) {
     std::cerr << "\nMixMaxRng::getState(): vector has wrong checksum!"
               << "\nInput vector is probably mispositioned now.\n";
     return false;
   }
   return true;
}

myuint_t MixMaxRng ::MOD_MULSPEC(myuint_t k)
{
   switch (N)
   {
      case 17:
          return 0;
          break;
      case 8:
          return 0;
          break;
      case 240:
          return fmodmulM61( 0, SPECIAL , (k) );
          break;
      default:
          std::cerr << "MIXMAX ERROR: " << "Disallowed value of parameter N\n";
          std::terminate();
          break;
   }
}

myuint_t MixMaxRng::MULWU (myuint_t k)
{
   return (( (k)<<(SPECIALMUL) & M61) ^ ( (k) >> (BITS-SPECIALMUL))  );
}

myuint_t MixMaxRng::iterate_raw_vec(myuint_t* Y, myuint_t sumtotOld)
{
   // operates with a raw vector, uses known sum of elements of Y
   int i;
            
   myuint_t  tempP, tempV;
   Y[0] = ( tempV = sumtotOld);
   myuint_t sumtot = Y[0], ovflow = 0; // will keep a running sum of all new elements
   tempP = 0;              // will keep a partial sum of all old elements
   for (i=1; (i<N); i++)
   {
     myuint_t tempPO = MULWU(tempP);
     tempP = modadd(tempP, Y[i]);
     tempV = MIXMAX_MOD_MERSENNE(tempV+tempP+tempPO); // new Y[i] = old Y[i] + old partial * m
     Y[i] = tempV;
     sumtot += tempV; if (sumtot < tempV) {ovflow++;}
   }
   return MIXMAX_MOD_MERSENNE(MIXMAX_MOD_MERSENNE(sumtot) + (ovflow <<3 ));
}
        
myuint_t MixMaxRng::get_next()
{
   int i;
   i=S.counter;
            
   if ((i<=(N-1)) )
   {
     S.counter++;
     return S.V[i];
   }
   else
   {
     S.sumtot = iterate_raw_vec(S.V.data(), S.sumtot);
     S.counter=2;
     return S.V[1];
   }
}

myuint_t MixMaxRng::precalc()
{
   int i;
   myuint_t temp;
   temp = 0;
   for (i=0; i < N; i++){
     temp = MIXMAX_MOD_MERSENNE(temp + S.V[i]);
   }
   S.sumtot = temp;
   return temp;
}
            
double MixMaxRng::get_next_float_packbits()
{
   myuint_t Z=get_next();
   return convert1double(Z);
}
        
void MixMaxRng::seed_vielbein(unsigned int index)
{
   int i;
   if (index<N)
   {
     for (i=0; i < N; i++){
       S.V[i] = 0;
     }
     S.V[index] = 1;
   }
   else
   {
     std::terminate();
   }
   S.counter = N;  // set the counter to N if iteration should happen right away
   S.sumtot = 1;
}
        
#define MIXMAX_SEED_WAS_ZERO 0xFF02

void MixMaxRng::seed_spbox(myuint_t seed)
{
   // a 64-bit LCG from Knuth line 26, in combination with a bit swap is used to seed

   const myuint_t MULT64=6364136223846793005ULL;
   int i;
            
   myuint_t sumtot=0,ovflow=0;
   if (seed == 0)
   {
     fprintf(stderr, " try seeding with nonzero seed next time!\n");
     exit(MIXMAX_SEED_WAS_ZERO);
   }
            
   myuint_t l = seed;
            
   for (i=0; i < N; i++){
     l*=MULT64; l = (l << 32) ^ (l>>32);
     S.V[i] = l & M61;
     sumtot += S.V[(i)]; if (sumtot < S.V[(i)]) {ovflow++;}
   }
   S.counter = N;  // set the counter to N if iteration should happen right away
   S.sumtot = MIXMAX_MOD_MERSENNE(MIXMAX_MOD_MERSENNE(sumtot) + (ovflow <<3 ));
}

#undef MIXMAX_SEED_WAS_ZERO

void MixMaxRng::seed_uniquestream( myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID )
{
   seed_vielbein(0);
   S.sumtot = apply_bigskip(S.V.data(), S.V.data(),  clusterID,  machineID,  runID,   streamID );
   S.counter = 1;
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
             
    Caution: never apply this to a derived vector, just choose some mother vector Vin, for example the unit vector by seed_vielbein(X,0),
    and use it in all your runs, just change runID to get completely nonoverlapping streams of random numbers on a different day.
             
    clusterID and machineID are provided for the benefit of large organizations who wish to ensure that a simulation
    which is running in parallel on a large number of  clusters and machines will have non-colliding source of random numbers.
             
    did i repeat it enough times? the non-collision guarantee is absolute, not probabilistic

   */

   const myuint_t skipMat17[128][17] =
   #include "CLHEP/Random/mixmax_skip_N17.icc"
   ;
            
   const myuint_t* skipMat[128];
   for (int i=0; i<128; i++) { skipMat[i] = skipMat17[i];}
            
   myID_t IDvec[4] = {streamID, runID, machineID, clusterID};
   int r,i,j,  IDindex;
   myID_t id;
   myuint_t Y[N], cum[N];
   myuint_t coeff;
   myuint_t* rowPtr;
   myuint_t sumtot=0;
            
   for (i=0; i<N; i++) { Y[i] = Vin[i]; sumtot = modadd( sumtot, Vin[i]); } ;
   for (IDindex=0; IDindex<4; IDindex++)
   { // go from lower order to higher order ID
     id=IDvec[IDindex];
     //printf("now doing ID at level %d, with ID = %d\n", IDindex, id);
     r = 0;
     while (id)
     {
       if (id & 1)
       {
         rowPtr = (myuint_t*)skipMat[r + IDindex*8*sizeof(myID_t)];
         for (i=0; i<N; i++){ cum[i] = 0; }
         for (j=0; j<N; j++)
         { // j is lag, enumerates terms of the poly
           // for zero lag Y is already given
           coeff = rowPtr[j]; // same coeff for all i
           for (i =0; i<N; i++){
             cum[i] =  fmodmulM61( cum[i], coeff ,  Y[i] ) ;
           }
           sumtot = iterate_raw_vec(Y, sumtot);
         }
         sumtot=0;
         for (i=0; i<N; i++){ Y[i] = cum[i]; sumtot = modadd( sumtot, cum[i]); } ;
       }
       id = (id >> 1); r++; // bring up the r-th bit in the ID
     }
   }
   sumtot=0;
   for (i=0; i<N; i++){ Vout[i] = Y[i]; sumtot = modadd( sumtot, Y[i]); }
   // returns sumtot, and copy the vector over to Vout
   return (sumtot) ;
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

myuint_t MixMaxRng::modadd(myuint_t foo, myuint_t bar)
{
#if defined(__x86_64__) && defined(__GNUC__) && (!defined(__ICC))
   //#warning Using assembler routine in modadd
   myuint_t out;
   /* Assembler trick suggested by Andrzej Görlich     */
   __asm__ ("addq %2, %0; "
            "btrq $61, %0; "
            "adcq $0, %0; "
            :"=r"(out)
            :"0"(foo), "r"(bar)
           );
   return out;
#else
   return MIXMAX_MOD_MERSENNE(foo+bar);
#endif
}

void MixMaxRng::print_state() const
{
   int j;
   std::cout << "mixmax state, file version 1.0\n";
   std::cout << "N=" << rng_get_N() << "; V[N]={";
   for (j=0; (j< (rng_get_N()-1) ); j++) {
     std::cout << S.V[j] << ", ";
   }
   std::cout << S.V[rng_get_N()-1];
   std::cout << "}; ";
   std::cout << "counter= " << S.counter;
   std::cout << "sumtot= " << S.sumtot << "\n";
}

MixMaxRng MixMaxRng::Branch()
{
   S.sumtot = iterate_raw_vec(S.V.data(), S.sumtot); S.counter = 1;
   MixMaxRng tmp=*this;
   tmp.BranchInplace(0); // daughter id
   return tmp;
}
    
void MixMaxRng::BranchInplace(int id)
{
   // Dont forget to iterate the mother, when branching the daughter, or else will have collisions!
   // a 64-bit LCG from Knuth line 26, is used to mangle a vector component
   constexpr myuint_t MULT64=6364136223846793005ULL;
   myuint_t tmp=S.V[id];
   S.V[1] *= MULT64; S.V[id] &= M61;
   S.sumtot = MIXMAX_MOD_MERSENNE( S.sumtot + S.V[id] - tmp + M61);
   S.sumtot = iterate_raw_vec(S.V.data(), S.sumtot);// printf("iterating!\n");
   S.counter = 1;
}

}  // namespace CLHEP
