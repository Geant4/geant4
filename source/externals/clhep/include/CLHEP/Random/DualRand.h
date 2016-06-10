// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                           Hep Random
//                        --- DualRand ---
//                        class header file
// -----------------------------------------------------------------------
//									
//  Canopy random number generator DualRand
//	Re-written as C++ routine for 32-bit ints MF 1/26/98
//
//  Exclusive or of a feedback shift register and integer congruence
//  random number generator.  The feedback shift register uses offsets
//  127 and 97.  The integer congruence generator uses a different
//  multiplier for each stream.  The multipliers are chosen to give
//  full period and maximum "potency" for modulo 2^32.  The period of
//  the combined random number generator is 2^159 - 2^32, and the
//  sequences are different for each stream (not just started in a
//  different place).
//
// =======================================================================
//  Canopy random number generator DualRand.
//  Doug Toussaint   5/25/88					           
//  Optimized by GMH 7/26/88					           
//  Optimized by GMH 7/26/88					           
//  Repaired  by GMH 12/1/88 to update modular congruence state            
//  Put into ranlib by GMH 6/23/89				           
//  Re-written as C++ routine for 32-bit ints MF 1/26/98	           
//  Re-written for CLHEP package	     KLS 6/04/98	           
//  Removed pow() from flat method for speed KLS 7/21/98	           
//  Ken Smith	   - Added conversion operators:  6th Aug 1998             
//  Mark Fischler    methods for distrib. instance save/restore 12/8/04    
//  Mark Fischler    methods for anonymous save/restore 12/27/04    
// Mark Fischler  - methods for vector save/restore 3/7/05    
// =======================================================================


#ifndef DualRand_h
#define DualRand_h

#include "CLHEP/Random/RandomEngine.h"

namespace CLHEP {

/**
 * @author
 * @ingroup random
 */
class DualRand: public HepRandomEngine {

public:

  DualRand();
  DualRand(long seed);
  DualRand(std::istream & is);
  DualRand(int rowIndex, int colIndex);
  virtual ~DualRand();

  // let the compiler generate the copy constructors
  //DualRand(const DualRand & p);
  //DualRand & operator=(const DualRand & p);

  double flat();
  // Returns a pseudo random number between 0 and 1 
  // (excluding the end points)

  void flatArray(const int size, double * vect);
  // Fills an array "vect" of specified size with flat random values.

  void setSeed(long seed, int);
  // Sets the state of the algorithm according to seed.

  void setSeeds(const long * seeds, int);
  // Sets the state of the algorithm according to the zero-terminated 
  // array of seeds.

  void saveStatus( const char filename[] = "DualRand.conf") const;
  // Saves on named file the current engine status.

  void restoreStatus( const char filename[] = "DualRand.conf" );
  // Reads from named file the last saved engine status and restores it.

  void showStatus() const;
  // Dumps the current engine status on the screen.

  operator float();      // flat value, without worrying about filling bits
  operator unsigned int();  // 32-bit flat value, quickest of all

  virtual std::ostream & put (std::ostream & os) const;
  virtual std::istream & get (std::istream & is);
  static  std::string beginTag ( );
  virtual std::istream & getState ( std::istream & is );
    
  std::string name() const;
  static std::string engineName() {return "DualRand";}

  std::vector<unsigned long> put () const;
  bool get (const std::vector<unsigned long> & v);
  bool getState (const std::vector<unsigned long> & v);
  
  static const unsigned int VECTOR_STATE_SIZE = 9;

private:

  // This generator is composed of two others combined:

  class Tausworthe {
  public:
    Tausworthe();
    Tausworthe(unsigned int seed);
    operator unsigned int();
    void put(std::ostream & os) const;
    void put(std::vector<unsigned long> & v) const;
    void get(std::istream & is);
    bool get(std::vector<unsigned long>::const_iterator & iv);
  private:
    int wordIndex;
    unsigned int words[4];
  }; // Tausworthe

  class IntegerCong {
  public:
    IntegerCong();
    IntegerCong(unsigned int seed, int streamNumber);
    operator unsigned int();
    void put(std::ostream & os) const;
    void put(std::vector<unsigned long> & v) const;
    void get(std::istream & is);
    bool get(std::vector<unsigned long>::const_iterator & iv);
  private:
    unsigned int state, multiplier, addend;
  }; // IntegerCong

  int numEngines;
  Tausworthe  tausworthe;
  IntegerCong integerCong;

}; // DualRand

}  // namespace CLHEP

#endif // DualRand_h
