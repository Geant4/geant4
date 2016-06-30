#ifndef CLHEP_ATOMIC_INT_H
#define CLHEP_ATOMIC_INT_H

// ======================================================================
//
// Use std::atomic when the compiler declares it uses the C++11 standard
//
// ======================================================================

#if (defined (G4MULTITHREADED))

  #if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 7)
    #include <atomic>
    #define CLHEP_ATOMIC_INT_TYPE std::atomic<int>
  #elif __clang__
    #if __has_feature(c_atomic)
      #include <atomic>
      #define CLHEP_ATOMIC_INT_TYPE std::atomic<int>
    #else
      #define CLHEP_ATOMIC_INT_TYPE int
    #endif
  #else
    #define CLHEP_ATOMIC_INT_TYPE int
  #endif

#else
  #define CLHEP_ATOMIC_INT_TYPE int
#endif

#endif
