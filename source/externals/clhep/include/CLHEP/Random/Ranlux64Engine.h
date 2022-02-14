// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                        --- Ranlux64Engine ---
//                          class header file
// -----------------------------------------------------------------------
// The algorithm for this random engine has been taken from the notes of 
// a double-precision ranlux implementation by Martin Luscher, dated 
// November 1997.
//
// Like the previous ranlux generator, this one also has "luxury" levels,
// determining how many pseudo-random numbers are discarded for every 
// twelve values used. Three levels are given, with the note that Luscher
// himself advocates only the highest two levels for this engine.
//  level 0  (p=109):              Throw away 109 values for every 12 used
//  level 1  (p=202):   (default)  Throw away 202 values for every 12 used
//  level 2  (p=397):              Throw away 397 values for every 12 used
//
// The initialization is carried out using a Multiplicative Congruential
// generator using formula constants of L'Ecuyer as described in "F.James,
// Comp. Phys. Comm. 60 (1990) 329-344".
// =======================================================================
// Ken Smith      - Created Initial draft: 14th Jul 1998
//                - Added conversion operators:  6th Aug 1998
// Mark Fischler
// 9/9/98	  - Added update() routine to allow computation of many at once
//		  - Replaced algorithm with jone exactly matching Luscher:
//			48-bits generated
//			skip n-12 instead of n numbers
//		  - Corrected protection agains overflow
// 12/8/04        - Methods for instance save/restore     
// 12/27/04       - methods for anonymous save/restore 12/27/04    
//
// =======================================================================

#ifndef Ranlux64Engine_h
#define Ranlux64Engine_h

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {

/**
 * @author
 * @ingroup random
 */
class Ranlux64Engine : public HepRandomEngine {

public:

  Ranlux64Engine( std::istream& is );
  Ranlux64Engine();
  Ranlux64Engine( long seed, int lxr = 1 );
  Ranlux64Engine( int rowIndex, int colIndex, int lxr );
  virtual ~Ranlux64Engine();
  // Constructors and destructor

  double flat();
  // It returns a pseudo random number between 0 and 1,
  // excluding the end points.

  void flatArray (const int size, double* vect);
  // Fills the array "vect" of specified size with flat random values.

  void setSeed(long seed, int lxr=1);
  // Sets the state of the algorithm according to seed.

  void setSeeds(const long * seeds, int lxr=1);
  // Sets the state of the algorithm according to the zero terminated
  // array of seeds.  Only the first seed is used.

  void saveStatus( const char filename[] = "Ranlux64.conf" ) const;
  // Saves in named file the current engine status.

  void restoreStatus( const char filename[] = "Ranlux64.conf" );
  // Reads from named file the last saved engine status and restores it.

  void showStatus() const;
  // Dumps the engine status on the screen.

  int getLuxury() const { return luxury; }
  // Gets the luxury level.

  virtual std::ostream & put (std::ostream & os) const;
  virtual std::istream & get (std::istream & is);
  static  std::string beginTag ( );
  virtual std::istream & getState ( std::istream & is );

  std::string name() const;
  static std::string engineName() {return "Ranlux64Engine";}

  std::vector<unsigned long> put () const;
  bool get (const std::vector<unsigned long> & v);
  bool getState (const std::vector<unsigned long> & v);
  
  static const unsigned int VECTOR_STATE_SIZE = 30;
  
private:

  void update();
  void advance(int dozens);

  int pDiscard;     // separate sequence by p-r = p-12 discarded elements
  int pDozens;      // pDiscard / 12;
  int endIters;     // pDiscard % 12;
  int luxury;

  int index;
  double randoms[12]; // randoms [i] is the x[n-i] of Luscher's note
  double carry;

}; // Ranlux64Engine

}  // namespace CLHEP

#endif // Ranlux64Engine_h
