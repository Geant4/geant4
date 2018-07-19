#ifndef CLHEP_USE_ATOMIC_GUARD_H
#define CLHEP_USE_ATOMIC_GUARD_H

// ======================================================================
//
// Use std::atomic when the compiler declares it uses the C++11 standard
//
// ======================================================================

#if defined (G4MULTITHREADED)

  #if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 7) || __clang__ || WIN32
    #include <atomic>
    #define CLHEP_USE_ATOMIC
  #endif

#endif

#endif
