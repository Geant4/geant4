#ifndef HEP_DEFS_H
#define HEP_DEFS_H

#ifdef WIN32
  //
  // Define DLL export macro for WIN32 systems for
  // importing/exporting external symbols to DLLs
  //
  #if defined G4LIB_BUILD_DLL
    #if defined CLHEP_EXPORT
      #define DLL_API __declspec( dllexport )
    #else
      #define DLL_API __declspec( dllimport )
    #endif
  #else
    #define DLL_API
  #endif
#else
  #define DLL_API
#endif

#endif /* HEP_DEFS_H */
