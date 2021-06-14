// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                        --- RanluxEngine ---
//                          class header file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).
//
// The algorithm for this random engine has been taken from the original
// implementation in FORTRAN by Fred James as part of the MATHLIB HEP
// library.
// The initialisation is carried out using a Multiplicative Congruential
// generator using formula constants of L'Ecuyer as described in "F.James,
// Comp. Phys. Comm. 60 (1990) 329-344".

// =======================================================================
// Adeyemi Adesanya - Created: 6th November 1995
// Gabriele Cosmo - Adapted & Revised: 22nd November 1995
// Adeyemi Adesanya - Added setSeeds() method: 2nd February 1996
// Gabriele Cosmo - Added flatArray() method: 8th February 1996
//                - Added methods for engine status: 19th November 1996
//                - Added default luxury value for setSeed()
//                  and setSeeds(): 21st July 1997
// J.Marraffino   - Added stream operators and related constructor.
//                  Added automatic seed selection from seed table and
//                  engine counter: 14th Feb 1998
// Ken Smith      - Added conversion operators:  6th Aug 1998
// Mark Fischler    Methods put, get for instance save/restore 12/8/04    
// Mark Fischler    methods for anonymous save/restore 12/27/04    
// =======================================================================

#ifndef RanluxEngine_h
#define RanluxEngine_h 1

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {

/**
 * @author
 * @ingroup random
 */
class RanluxEngine : public HepRandomEngine {

public:

  RanluxEngine( std::istream& is );
  RanluxEngine();
  RanluxEngine( long seed, int lxr = 3 );
  RanluxEngine( int rowIndex, int colIndex, int lxr );
  virtual ~RanluxEngine();
  // Constructors and destructor

// Luxury level is set in the same way as the original FORTRAN routine.
//  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
//           and Zaman, very long period, but fails many tests.
//  level 1  (p=48): considerable improvement in quality over level 0,
//           now passes the gap test, but still fails spectral test.
//  level 2  (p=97): passes all known tests, but theoretically still
//           defective.
//  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
//           correlations have very small chance of being observed.
//  level 4  (p=389): highest possible luxury, all 24 bits chaotic.

  double flat();
  // It returns a pseudo random number between 0 and 1,
  // excluding the end points.

  void flatArray (const int size, double* vect);
  // Fills the array "vect" of specified size with flat random values.

  void setSeed(long seed, int lxr=3);
  // Sets the state of the algorithm according to seed.

  void setSeeds(const long * seeds, int lxr=3);
  // Sets the state of the algorithm according to the zero terminated
  // array of seeds. Only the first seed is used.

  void saveStatus( const char filename[] = "Ranlux.conf" ) const;
  // Saves on file Ranlux.conf the current engine status.

  void restoreStatus( const char filename[] = "Ranlux.conf" );
  // Reads from file Ranlux.conf the last saved engine status
  // and restores it.

  void showStatus() const;
  // Dumps the engine status on the screen.

  int getLuxury() const { return luxury; }
  // Gets the luxury level.

  operator double();       // Returns same as flat()
  operator float();        // less precise flat, faster if possible
  operator unsigned int(); // 32-bit flat, but slower than double or float

  virtual std::ostream & put (std::ostream & os) const;
  virtual std::istream & get (std::istream & is);
  static  std::string beginTag ( );
  virtual std::istream & getState ( std::istream & is );

  std::string name() const;
  static std::string engineName() {return "RanluxEngine";}

  std::vector<unsigned long> put () const;
  bool get (const std::vector<unsigned long> & v);
  bool getState (const std::vector<unsigned long> & v);
  
  static const unsigned int VECTOR_STATE_SIZE = 31;
  
private:

  int nskip, luxury;
  float float_seed_table[24];
  int i_lag,j_lag;  
  float carry;
  int count24;
  static const int int_modulus = 0x1000000;
};

}  // namespace CLHEP

#endif
