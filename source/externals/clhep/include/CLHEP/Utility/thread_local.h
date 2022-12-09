#ifndef CLHEP_THREAD_LOCAL_H
#define CLHEP_THREAD_LOCAL_H

// ======================================================================
//
// Use thread_local when the compiler declares it uses the C++11 standard
//
// ======================================================================

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 7) || __clang__ || defined(WIN32) || defined(__MINGW32__)
  #define CLHEP_THREAD_LOCAL thread_local
#else
  #define CLHEP_THREAD_LOCAL
#endif

#endif
