// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                          HEP Random
//                       --- MixMaxRng ---
//                       class header file
// -----------------------------------------------------------------------
//
// This file interfaces the PseudoRandom Number Generator 
// proposed by:
// N.Z. Akopov, G.K.Saviddy & N.G. Ter-Arutyunian
//   "Matrix Generator of Pseudorandom Numbers",
//   J.Compt.Phy. 97, 573 (1991)
//   Preprint: EPI-867(18)-86, Yerevan June 1986.
// G. Savvidy & N. Savvidy
//   "On the Monte Carlo Simulation of Physical Systems",
//  J.Comput.Phys. 97 (1991) 566

// =======================================================================
// Implementation by Konstantin Savvidy - 2004-2015
//  Release 0.99 and later: released under the LGPL license version 3.0
// =======================================================================
// CLHEP interface implemented by 
//     J. Apostolakis, G. Cosmo & K. Savvidy - Created: 6th July 2015
//     CLHEP interface released under the LGPL license version 3.0
// =======================================================================

#ifndef MixMaxRng_h
#define MixMaxRng_h 1

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/mixmax.h"

namespace CLHEP {

/**
 * @author  K. Savvidy
 * @ingroup random
 */
class MixMaxRng: public HepRandomEngine {

public:

  MixMaxRng(std::istream& is);
  MixMaxRng();
  MixMaxRng(long seed);
  MixMaxRng(int rowIndex, int colIndex);
  virtual ~MixMaxRng();
  // Constructor and destructor.

  MixMaxRng(const MixMaxRng& rng);
  MixMaxRng& operator=(const MixMaxRng& rng);
  // Copy constructor and assignment operator.

  double flat();
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

  operator unsigned int();
  // 32-bit flat

  virtual std::ostream & put (std::ostream & os) const;
  virtual std::istream & get (std::istream & is);
  static  std::string beginTag ( );
  virtual std::istream & getState ( std::istream & is );

  std::string name() const;
  static std::string engineName() {return "MixMaxRng";}

  std::vector<unsigned long> put () const;
  bool get (const std::vector<unsigned long> & v);
  bool getState (const std::vector<unsigned long> & v);
  
  static const unsigned int VECTOR_STATE_SIZE = 2*N+4; // 2N+4 for MIXMAX
  
private:

  // Pointer to the current status of the generator.
  rng_state_st* fRngState;
};

}  // namespace CLHEP

#endif
