// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                           HEP Random
//                      --- RanshiEngine ---
//                       class header file
// -----------------------------------------------------------------------
//
//
// The algorithm for this random engine was taken from "F.Gutbrod, Comp.
// Phys. Comm. 87 (1995) 291-306".
//
// The algorithm can be imagined as a physical system as follows: Imagine
// 512 "black balls" each with their own unique spin, and positions char- 
// acterized by disrete angles, where the spin is a 32-bit unsigned integer.
// A "red ball" collides based upon the angle determined by the last 8 bits
// of its spin, and the spin of the colliding ball is taken as the output
// random number. The spin of the colliding ball is replaced then with the
// left circular shift of the black ball's spin XOR'd with the red ball's
// spin. The black ball's old spin becomes the red ball's.
//
// To avoid the traps presented, two measures are taken: first, the red 
// ball will oscillate between hitting the lower half of the buffer on one
// turn and the upper half on another; second, the red ball's spin is 
// incremented by a counter of the number of random numbers produced.
//
// The result is scaled to a double precision floating point number to which
// is added another random double further scaled 2^(53-32) places to the
// right in order to ensure that the remaining bits of the result are not 
// left empty due to the mere 32 bits representation used internally.

// =======================================================================
// Ken Smith      - Created: 9th June 1998
//                - Removed pow() from flat method: 21st Jul 1998
//                - Added conversion operators:  6th Aug 1998
// Mark Fischler    Methods put, get for instance save/restore 12/8/04    
// Mark Fischler    methods for anonymous save/restore 12/27/04    
// =======================================================================

#ifndef HepRanshiEngine_h
#define HepRanshiEngine_h

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {

/**
 * @author
 * @ingroup random
 */
class RanshiEngine: public HepRandomEngine {
 
public:

    RanshiEngine();
    RanshiEngine(std::istream &is);
    RanshiEngine(long seed);
    RanshiEngine(int rowIndex, int colIndex);
    virtual ~RanshiEngine();
    // Constructors and destructor

    double flat();
    // Returns a pseudo random number between 0 and 1

    void flatArray(const int size, double* vect);
    // Fills the array "vect" of specified size with flat random values

    void setSeed(long seed, int);
    // Sets the state of the algorithm according to seed.

    void setSeeds(const long* seeds, int);
    // Sets the state of the algorithm according to the zero-terminated
    // array of seeds. 

    void saveStatus(const char filename[] = "RanshiEngine.conf") const;
    // Saves on named file the current engine status

    void restoreStatus(const char filename[] = "RanshiEngine.conf");
    // Reads from named file the last saved engine status
    // and restores it.

    void showStatus() const;
    // Dumps the engine status on the screen

    operator float();      // flat value, without worrying about filling bits
    operator unsigned int();  // 32-bit flat value, quickest of all

   virtual std::ostream & put (std::ostream & os) const;
   virtual std::istream & get (std::istream & is);
  static  std::string beginTag ( );
  virtual std::istream & getState ( std::istream & is );

   std::string name() const;
   static std::string engineName() {return "RanshiEngine";}

  std::vector<unsigned long> put () const;
  bool get (const std::vector<unsigned long> & v);
  bool getState (const std::vector<unsigned long> & v);
  
private:
    enum {numBuff = 512};

    unsigned int halfBuff, numFlats;
    unsigned int buffer[numBuff];
    unsigned int redSpin;

    static const unsigned int VECTOR_STATE_SIZE = numBuff + 4;
  
}; // RanshiEngine

}  // namespace CLHEP

#endif // HepRanshiEngine_h
