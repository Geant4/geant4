//
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                          HEP Random
//                       --- MixMaxRng ---
//                       class header file
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

#ifndef MixMaxRng_h
#define MixMaxRng_h 1

#include <array>
#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {
    
/**
  * @author  K.Savvidy
  * @ingroup random
  */

using myID_t = std::uint32_t;
using myuint_t = unsigned long long int;

class MixMaxRng: public HepRandomEngine {

  static const int N = 17;

public:

  MixMaxRng(std::istream& is);
  MixMaxRng();
  MixMaxRng(long seed);
  ~MixMaxRng();
  // Constructors and destructor.

  MixMaxRng(const MixMaxRng& rng);
  MixMaxRng& operator=(const MixMaxRng& rng);
  // Copy constructor and assignment operator.

  double flat() { return (S.counter<=(N-1)) ? generate(S.counter):iterate(); }
  // Returns a pseudo random number between 0 and 1
  // (excluding the zero: in (0,1] )
  // smallest number which it will give is approximately 10^-19

  void flatArray (const int size, double* vect);
  // Fills the array "vect" of specified size with flat random values.

  void setSeed(long seed, int dum=0);
  // Sets the state of the algorithm according to seed.

  void setSeeds(const long * seeds, int seedNum=0);
  // Sets the initial state of the engine according to the array of between one and four 32-bit seeds.
  // If the size of long is greater on the platform, only the lower 32-bits are used.
  // Streams created from seeds differing by at least one bit somewhere are guaranteed absolutely
  // to be independent and non-colliding for at least the next 10^100 random numbers

  void saveStatus( const char filename[] = "MixMaxRngState.conf" ) const;
  // Saves the the current engine state in the file given, by default MixMaxRngState.conf

  void restoreStatus( const char filename[] = "MixMaxRngState.conf" );
  // Reads a valid engine state from a given file, by default MixMaxRngState.conf
  // and restores it.

  void showStatus() const;
  // Dumps the engine status on the screen.

  operator double();
  // Returns same as flat()
  operator float();
  // less precise flat, faster if possible
  operator unsigned int();
  // 32-bit flat

  virtual std::ostream & put (std::ostream & os) const;
  virtual std::istream & get (std::istream & is);
  static  std::string beginTag ( );
  virtual std::istream & getState ( std::istream & is );

  std::string name() const { return "MixMaxRng"; }
  static std::string engineName();

  std::vector<unsigned long> put () const;
  bool get (const std::vector<unsigned long> & v);
  bool getState (const std::vector<unsigned long> & v);

private:

  static constexpr long long int SPECIAL   = ((N==17)? 0 : ((N==240)? 487013230256099140ULL:0) ); // etc...
  static constexpr long long int SPECIALMUL= ((N==17)? 36: ((N==240)? 51                   :53) ); // etc...
  // Note the potential for confusion...
  static constexpr int BITS=61;
  static constexpr myuint_t M61=2305843009213693951ULL;
  static constexpr double INV_M61=0.43368086899420177360298E-18;
  static constexpr unsigned int VECTOR_STATE_SIZE = 2*N+4; // 2N+4 for MIXMAX

  #define MIXMAX_MOD_MERSENNE(k) ((((k)) & M61) + (((k)) >> BITS) )

  static constexpr int rng_get_N();
  static constexpr long long int rng_get_SPECIAL();
  static constexpr int rng_get_SPECIALMUL();
  void seed_uniquestream( myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID );
  void seed_spbox(myuint_t);
  void print_state() const;
  myuint_t precalc();
  myuint_t get_next() ;
  inline double get_next_float() { return get_next_float_packbits(); }
  // Returns a random double with all 52 bits random, in the range (0,1]
  
  MixMaxRng Branch();
  void BranchInplace(int id);

  MixMaxRng(myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID );	   // Constructor with four 32-bit seeds
  inline void seed64(myuint_t seedval) { seed_uniquestream( 0, 0, (myID_t)(seedval>>32), (myID_t)seedval ); } // seed with one 64-bit seed

  double generate(int i);
  double iterate();

  double get_next_float_packbits();
#if defined __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif
  inline double convert1double(myuint_t u)
  {
    const double one = 1;
    const myuint_t onemask = *(myuint_t*)&one;
    myuint_t tmp = (u>>9) | onemask; // bits between 52 and 62 dont affect the result!
    double d = *(double*)&tmp;
    return d-1.0;
  }
#if defined __GNUC__
#pragma GCC diagnostic pop
#endif
  myuint_t MOD_MULSPEC(myuint_t k);
  myuint_t MULWU(myuint_t k);
  void seed_vielbein( unsigned int i); // seeds with the i-th unit vector, i = 0..N-1,  for testing only
  myuint_t iterate_raw_vec(myuint_t* Y, myuint_t sumtotOld);
  myuint_t apply_bigskip(myuint_t* Vout, myuint_t* Vin, myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID );
  myuint_t modadd(myuint_t foo, myuint_t bar);
#if defined(__x86_64__)
  myuint_t mod128(__uint128_t s);
  myuint_t fmodmulM61(myuint_t cum, myuint_t a, myuint_t b);
#else // on all other platforms, including 32-bit linux, PPC and PPC64, ARM and all Windows
  myuint_t fmodmulM61(myuint_t cum, myuint_t s, myuint_t a);
#endif

private:

  struct rng_state_st
  {
      std::array<myuint_t, N> V;
      myuint_t sumtot;
      int counter;
  };

  typedef struct rng_state_st rng_state_t;     // struct alias
  rng_state_t S;
};

}  // namespace CLHEP

#endif
