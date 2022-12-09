//
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                       --- RanluxppEngine ---
//                         class header file
// -----------------------------------------------------------------------
// Implementation of the RANLUX++ generator
//
// RANLUX++ is an LCG equivalent of RANLUX using 576 bit numbers.
//
// References:
// A. Sibidanov
//   A revision of the subtract-with-borrow random numbergenerators
//   Computer Physics Communications, 221(2017), 299-303
//
// J. Hahnfeld, L. Moneta
//   A Portable Implementation of RANLUX++
//   vCHEP2021

#ifndef RanluxppEngine_h
#define RanluxppEngine_h

#include "CLHEP/Random/RandomEngine.h"

#include <cstdint>

namespace CLHEP {

/**
 * @author Jonas Hahnfeld
 * @ingroup random
 */
class RanluxppEngine final : public HepRandomEngine {

public:
  RanluxppEngine();
  RanluxppEngine(long seed);
  RanluxppEngine(std::istream &is);
  virtual ~RanluxppEngine();
  // Constructors and destructor

  double flat() override;
  // It returns a pseudo random number between 0 and 1,
  // excluding the end points.

  void flatArray(const int size, double *vect) override;
  // Fills the array "vect" of specified size with flat random values.

  void setSeed(long seed, int dummy = 0) override;
  // Sets the state of the algorithm according to seed.

  void setSeeds(const long *seeds, int dummy = 0) override;
  // Sets the state of the algorithm according to the zero terminated
  // array of seeds.  Only the first seed is used.

  void skip(uint64_t n);
  // Skip `n` random numbers without generating them.

  void saveStatus(const char filename[] = "Ranluxpp.conf") const override;
  // Saves in named file the current engine status.

  void restoreStatus(const char filename[] = "Ranluxpp.conf") override;
  // Reads from named file the last saved engine status and restores it.

  void showStatus() const override;
  // Dumps the engine status on the screen.

  std::string name() const override;

  // Optional methods to serialize the engine's state into vectors and streams.
  static std::string engineName();
  static std::string beginTag();

  std::ostream &put(std::ostream &os) const override;
  std::istream &get(std::istream &is) override;

  std::istream &getState(std::istream &is) override;

  std::vector<unsigned long> put() const override;
  bool get(const std::vector<unsigned long> &v) override;
  bool getState(const std::vector<unsigned long> &v) override;

  // Save and restore to/from streams
  operator double() override { return flat(); }
  operator float() override { return float(flat()); }
  operator unsigned int() override { return (unsigned int)nextRandomBits(); }

  // 1 value for the engine ID, 2 * 9 values for the state, and 2 more values
  // for the carry bit and the position.
  static const unsigned int VECTOR_STATE_SIZE = 21;

private:
  void advance();
  uint64_t nextRandomBits();

  uint64_t fState[9]; ///< RANLUX state of the generator
  unsigned fCarry;    ///< Carry bit of the RANLUX state
  int fPosition = 0;  ///< Current position in bits

}; // RanluxppEngine

} // namespace CLHEP

#endif
