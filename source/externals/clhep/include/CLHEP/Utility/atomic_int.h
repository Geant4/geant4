#ifndef CLHEP_ATOMIC_INT_H
#define CLHEP_ATOMIC_INT_H

// ======================================================================
//
// Use std::atomic when the compiler declares it uses the C++11 standard
//
// ======================================================================

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 7) || __clang__ || WIN32
  #include <atomic>
  #define CLHEP_ATOMIC_INT_TYPE std::atomic<int>
#else
  #define CLHEP_ATOMIC_INT_TYPE int
#endif

#endif
