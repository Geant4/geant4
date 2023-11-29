//
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                       --- RanluxppEngine ---
//                     class implementation file
// -----------------------------------------------------------------------
// Implementation of the RANLUX++ generator
//
// RANLUX++ is an LCG equivalent of RANLUX using 576 bit numbers.
//
// The idea of the generator (such as the initialization method) and the algorithm
// for the modulo operation are described in
// A. Sibidanov, *A revision of the subtract-with-borrow random numbergenerators*,
// *Computer Physics Communications*, 221(2017), 299-303,
// preprint https://arxiv.org/pdf/1705.03123.pdf
//
// The code is loosely based on the Assembly implementation by A. Sibidanov
// available at https://github.com/sibidanov/ranluxpp/.
//
// Compared to the original generator, this implementation contains a fix to ensure
// that the modulo operation of the LCG always returns the smallest value congruent
// to the modulus (based on notes by M. Lüscher). Also, the generator converts the
// LCG state back to RANLUX numbers (implementation based on notes by M. Lüscher).
// This avoids a bias in the generated numbers because the upper bits of the LCG
// state, that is smaller than the modulus \f$ m = 2^{576} - 2^{240} + 1 \f$ (not
// a power of 2!), have a higher probability of being 0 than 1. And finally, this
// implementation draws 48 random bits for each generated floating point number
// (instead of 52 bits as in the original generator) to maintain the theoretical
// properties from understanding the original transition function of RANLUX as a
// chaotic dynamical system.
//
// These modifications and the portable implementation in general are described in
// J. Hahnfeld, L. Moneta, *A Portable Implementation of RANLUX++*, vCHEP2021
// preprint https://arxiv.org/pdf/2106.02504.pdf

#include "CLHEP/Random/RanluxppEngine.h"

#include "CLHEP/Random/engineIDulong.h"
#include "CLHEP/Utility/atomic_int.h"

#include "CLHEP/Random/ranluxpp/mulmod.h"
#include "CLHEP/Random/ranluxpp/ranlux_lcg.h"

#include <cassert>
#include <fstream>
#include <ios>
#include <cstdint>

namespace CLHEP {

namespace {
// Number of instances with automatic seed selection.
CLHEP_ATOMIC_INT_TYPE numberOfEngines(0);

const uint64_t kA_2048[] = {
    0xed7faa90747aaad9, 0x4cec2c78af55c101, 0xe64dcb31c48228ec,
    0x6d8a15a13bee7cb0, 0x20b2ca60cb78c509, 0x256c3d3c662ea36c,
    0xff74e54107684ed2, 0x492edfcc0cc8e753, 0xb48c187cf5b22097,
};
} // namespace

RanluxppEngine::RanluxppEngine() : HepRandomEngine() {
  int numEngines = ++numberOfEngines;
  setSeed(numEngines);
}

RanluxppEngine::RanluxppEngine(long seed) : HepRandomEngine() {
  theSeed = seed;
  setSeed(seed);
}

RanluxppEngine::RanluxppEngine(std::istream &is) : HepRandomEngine() {
  get(is);
}

RanluxppEngine::~RanluxppEngine() = default;

static constexpr int kMaxPos = 9 * 64;
static constexpr int kBits = 48;

void RanluxppEngine::advance() {
  uint64_t lcg[9];
  to_lcg(fState, fCarry, lcg);
  mulmod(kA_2048, lcg);
  to_ranlux(lcg, fState, fCarry);
  fPosition = 0;
}

uint64_t RanluxppEngine::nextRandomBits() {
  if (fPosition + kBits > kMaxPos) {
    advance();
  }

  int idx = fPosition / 64;
  int offset = fPosition % 64;
  int numBits = 64 - offset;

  uint64_t bits = fState[idx] >> offset;
  if (numBits < kBits) {
    bits |= fState[idx + 1] << numBits;
  }
  bits &= ((uint64_t(1) << kBits) - 1);

  fPosition += kBits;
  assert(fPosition <= kMaxPos && "position out of range!");

  return bits;
}

double RanluxppEngine::flat() {
  // RandomEngine wants a "double random values ranging between ]0,1[", so
  // exclude all zero bits.
  uint64_t random;
  do {
    random = nextRandomBits();
  } while (random == 0);

  static constexpr double div = 1.0 / (uint64_t(1) << kBits);
  return random * div;
}

void RanluxppEngine::flatArray(const int size, double *vect) {
  for (int i = 0; i < size; i++) {
    vect[i] = flat();
  }
}

void RanluxppEngine::setSeed(long seed, int) {
  theSeed = seed;

  uint64_t lcg[9];
  lcg[0] = 1;
  for (int i = 1; i < 9; i++) {
    lcg[i] = 0;
  }

  uint64_t a_seed[9];
  // Skip 2 ** 96 states.
  powermod(kA_2048, a_seed, uint64_t(1) << 48);
  powermod(a_seed, a_seed, uint64_t(1) << 48);
  // Skip more states according to seed.
  powermod(a_seed, a_seed, seed);
  mulmod(a_seed, lcg);

  to_ranlux(lcg, fState, fCarry);
  fPosition = 0;
}

void RanluxppEngine::setSeeds(const long *seeds, int) {
  theSeeds = seeds;
  setSeed(*seeds, 0);
}

void RanluxppEngine::skip(uint64_t n) {
  int left = (kMaxPos - fPosition) / kBits;
  assert(left >= 0 && "position was out of range!");
  if (n < (uint64_t)left) {
    // Just skip the next few entries in the currently available bits.
    fPosition += n * kBits;
    return;
  }

  n -= left;
  // Need to advance and possibly skip over blocks.
  int nPerState = kMaxPos / kBits;
  int skip = int(n / nPerState);

  uint64_t a_skip[9];
  powermod(kA_2048, a_skip, skip + 1);

  uint64_t lcg[9];
  to_lcg(fState, fCarry, lcg);
  mulmod(a_skip, lcg);
  to_ranlux(lcg, fState, fCarry);

  // Potentially skip numbers in the freshly generated block.
  int remaining = int(n - skip * nPerState);
  assert(remaining >= 0 && "should not end up at a negative position!");
  fPosition = remaining * kBits;
  assert(fPosition <= kMaxPos && "position out of range!");
}

void RanluxppEngine::saveStatus(const char filename[]) const {
  std::ofstream os(filename);
  put(os);
  os.close();
}

void RanluxppEngine::restoreStatus(const char filename[]) {
  std::ifstream is(filename);
  get(is);
  is.close();
}

void RanluxppEngine::showStatus() const {
  std::cout
      << "--------------------- RanluxppEngine status --------------------"
      << std::endl;
  std::cout << " fState[] = {";
  std::cout << std::hex << std::setfill('0');
  for (int i = 0; i < 9; i++) {
    if (i % 3 == 0) {
      std::cout << std::endl << "     ";
    } else {
      std::cout << " ";
    }
    std::cout << "0x" << std::setw(16) << fState[i] << ",";
  }
  std::cout << std::endl << " }" << std::endl;
  std::cout << std::dec;
  std::cout << " fCarry = " << fCarry << ", fPosition = " << fPosition
            << std::endl;
  std::cout
      << "----------------------------------------------------------------"
      << std::endl;
}

std::string RanluxppEngine::name() const { return engineName(); }

std::string RanluxppEngine::engineName() { return "RanluxppEngine"; }

std::string RanluxppEngine::beginTag() { return "RanluxppEngine-begin"; }

std::ostream &RanluxppEngine::put(std::ostream &os) const {
  os << beginTag() << "\n";
  const std::vector<unsigned long> state = put();
  for (unsigned long v : state) {
    os << v << "\n";
  }
  return os;
}

std::istream &RanluxppEngine::get(std::istream &is) {
  std::string tag;
  is >> tag;
  if (tag != beginTag()) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "No RanluxppEngine found at current position\n";
    return is;
  }
  return getState(is);
}

std::istream &RanluxppEngine::getState(std::istream &is) {
  std::vector<unsigned long> state;
  state.reserve(VECTOR_STATE_SIZE);
  for (unsigned int i = 0; i < VECTOR_STATE_SIZE; i++) {
    unsigned long v;
    is >> v;
    state.push_back(v);
  }

  getState(state);
  return is;
}

std::vector<unsigned long> RanluxppEngine::put() const {
  std::vector<unsigned long> v;
  v.reserve(VECTOR_STATE_SIZE);
  v.push_back(engineIDulong<RanluxppEngine>());

  // unsigned long is only guaranteed to be 32 bit wide, so chop up the 64 bit
  // values in fState.
  for (int i = 0; i < 9; i++) {
    unsigned long lower = static_cast<uint32_t>(fState[i]);
    v.push_back(lower);
    unsigned long upper = static_cast<uint32_t>(fState[i] >> 32);
    v.push_back(upper);
  }

  v.push_back(fCarry);
  v.push_back(fPosition);
  return v;
}

bool RanluxppEngine::get(const std::vector<unsigned long> &v) {
  if (v[0] != engineIDulong<RanluxppEngine>()) {
    std::cerr << "RanluxppEngine::get(): "
              << "vector has wrong ID word - state unchanged" << std::endl;
    return false;
  }
  return getState(v);
}

bool RanluxppEngine::getState(const std::vector<unsigned long> &v) {
  if (v.size() != VECTOR_STATE_SIZE) {
    std::cerr << "RanluxppEngine::getState(): "
              << "vector has wrong length - state unchanged" << std::endl;
    return false;
  }

  // Assemble the state vector (see RanluxppEngine::put).
  for (int i = 0; i < 9; i++) {
    uint64_t lower = v[2 * i + 1];
    uint64_t upper = v[2 * i + 2];
    fState[i] = (upper << 32) + lower;
  }
  fCarry = (unsigned int)v[19];
  fPosition = (int)v[20];

  return true;
}

} // namespace CLHEP
