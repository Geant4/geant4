#ifndef G3TOG4_DEFS_H
#define G3TOG4_DEFS_H

#ifdef WIN32
  //
  // Define DLL export macro for WIN32 systems for
  // importing/exporting external symbols to DLLs
  //
  #if defined G4LIB_BUILD_DLL
    #if defined G3TOG4_EXPORT
      #define G3G4DLL_API __declspec( dllexport )
    #else
      #define G3G4DLL_API __declspec( dllimport )
    #endif
  #else
    #define G3G4DLL_API
  #endif
#else
  #define G3G4DLL_API
#endif

#endif /* G3TOG4_DEFS_H */
