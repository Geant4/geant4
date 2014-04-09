// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                        --- HepJamesRandom ---
//                          class header file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).
//
// HepJamesRandom implements the algorithm by Marsaglia-Zaman RANMAR
// described in "F.James, Comp. Phys. Comm. 60 (1990) 329" and implemented
// in FORTRAN77 as part of the MATHLIB HEP library for pseudo-random
// numbers generation.
// This is the default random engine invoked by each distribution unless
// the user sets a different one.

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
//                - Minor corrections: 31st October 1996
//                - Added methods for engine status: 19th November 1996
//                - setSeed(), setSeeds() now have default dummy argument
//                  set to zero: 11th July 1997
// J.Marraffino   - Added stream operators and related constructor.
//                  Added automatic seed selection from seed table and
//                  engine counter: 16th Feb 1998
// Ken Smith      - Added conversion operators:  6th Aug 1998
// V. Innocente   - changed pointers to indices     3 may 2000
// Mark Fischler  - Methods for distrib. instance save/restore 12/8/04    
//  Mark Fischler    methods for anonymous save/restore 12/27/04    
// =======================================================================

#ifndef HepJamesRandom_h
#define HepJamesRandom_h 1

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {

/**
 * @author
 * @ingroup random
 */
class HepJamesRandom: public HepRandomEngine {

public:

  HepJamesRandom(std::istream& is);
  HepJamesRandom();
  HepJamesRandom(long seed);
  HepJamesRandom(int rowIndex, int colIndex);
  virtual ~HepJamesRandom();
  // Constructor and destructor.

  double flat();
  // Returns a pseudo random number between 0 and 1 
  // (excluding the end points)

  void flatArray (const int size, double* vect);
  // Fills the array "vect" of specified size with flat random values.

  void setSeed(long seed, int dum=0);
  // Sets the state of the algorithm according to seed.

  void setSeeds(const long * seeds, int dum=0);
  // Sets the state of the algorithm according to the zero terminated
  // array of seeds. Only the first seed is used.

  void saveStatus( const char filename[] = "JamesRand.conf" ) const;
  // Saves on file JamesRand.conf the current engine status.

  void restoreStatus( const char filename[] = "JamesRand.conf" );
  // Reads from file JamesRand.conf the last saved engine status
  // and restores it.

  void showStatus() const;
  // Dumps the engine status on the screen.

  operator unsigned int();
  // 32-bit flat, but slower than double or float.

  virtual std::ostream & put (std::ostream & os) const;
  virtual std::istream & get (std::istream & is);
  static  std::string beginTag ( );
  virtual std::istream & getState ( std::istream & is );

  std::string name() const;
  static std::string engineName() {return "HepJamesRandom";}

  std::vector<unsigned long> put () const;
  bool get (const std::vector<unsigned long> & v);
  bool getState (const std::vector<unsigned long> & v);
  
  static const unsigned int VECTOR_STATE_SIZE = 202;
  
private:

  // Members defining the current status of the generator.
  double u[97];
  double c, cd, cm;
  int i97, j97;
};

}  // namespace CLHEP

#endif
