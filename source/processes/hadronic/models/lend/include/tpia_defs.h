#ifndef tpia_defs_h_included
#define tpia_defs_h_included

#ifdef WIN32
  //
  // Define DLL export macro for WIN32 systems for
  // importing/exporting external symbols to DLLs
  //
  #if defined G4LIB_BUILD_DLL
    #if defined G4PROCESSES_EXPORT
      #define DLL_LEND __declspec( dllexport )
    #else
      #define DLL_LEND __declspec( dllimport )
    #endif
  #else
    #define DLL_LEND
  #endif
#else
  #define DLL_LEND
#endif

#endif          /* End of tpia_defs_h_included. */
