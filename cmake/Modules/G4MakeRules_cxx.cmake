# - Geant4MakeRules_cxx
# Sets the default make rules for a CXX build, specifically the
# initialization of the compiler flags on a platform and compiler
# dependent basis
#
#
#-----------------------------------------------------------------------
# function __configure_tls_models()
#          Set available thread local storage models. Valid for GNU,
#          Clang and Intel compilers. Adds an additional "auto"
#          dummy model to indicate that no flag should be added.
#
function(__configure_tls_models)
  # available models, default first
  set(_TLSMODELS initial-exec local-exec global-dynamic local-dynamic)
  foreach(_s ${_TLSMODELS})
    set(${_s}_TLSMODEL_FLAGS "-ftls-model=${_s}" PARENT_SCOPE)
  endforeach()

  list(APPEND _TLSMODELS auto)
  set(auto_TLSMODEL_FLAGS "" PARENT_SCOPE)

  set(TLSMODEL_IS_AVAILABLE ${_TLSMODELS} PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------
# DEFAULT FLAG SETTING
#-----------------------------------------------------------------------
# GNU C++ or Clang/AppleClang Compiler on all(?) platforms
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  # Warnings
  set(CMAKE_CXX_FLAGS_INIT "-W -Wall -pedantic -Wno-non-virtual-dtor -Wno-long-long -Wwrite-strings -Wpointer-arith -Woverloaded-virtual -Wno-variadic-macros -Wshadow")
  # Use pipes rather than temp files
  set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -pipe")
  # Remove superfluous "unused argument" and "GL deprecation" "warnings" from Clang
  if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -Qunused-arguments -DGL_SILENCE_DEPRECATION")
  endif()

  # Additional per-mode flags
  # Assist auto-vectorization in Release
  set(CMAKE_CXX_FLAGS_RELEASE_INIT "${CMAKE_CXX_FLAGS_RELEASE_INIT} -fno-trapping-math -ftree-vectorize -fno-math-errno")

  # Use Debug-safe optimization
  set(CMAKE_CXX_FLAGS_DEBUG_INIT "${CMAKE_CXX_FLAGS_DEBUG_INIT} -Og")

  # Extra Geant4 modes
  # - Debug_FPE: Core Debug mode plus FPE)
  set(CMAKE_CXX_FLAGS_DEBUG_FPE_INIT "${CMAKE_CXX_FLAGS_DEBUG_INIT} -DG4FPE_DEBUG")

  # - TestRelease
  set(CMAKE_CXX_FLAGS_TESTRELEASE_INIT "-g -DG4DEBUG_VERBOSE -DG4FPE_DEBUG")

  # - Maintainer
  set(CMAKE_CXX_FLAGS_MAINTAINER_INIT "-g")

  # - Multithreading
  __configure_tls_models()
  set(GEANT4_MULTITHREADED_CXX_FLAGS "-pthread")
endif()

#-----------------------------------------------------------------------
# MSVC - all (?) versions
if(MSVC)
  # Hmm, WIN32-VC.gmk uses dashes, but cmake uses slashes, latter probably
  # best for native build.
  set(CMAKE_CXX_FLAGS_INIT "-GR -EHsc -Zm200 -nologo -D_CONSOLE -D_WIN32 -DWIN32 -DOS -DXPNET -D_CRT_SECURE_NO_DEPRECATE")
  set(CMAKE_CXX_FLAGS_DEBUG_INIT "-MDd -Od -Zi")
  set(CMAKE_CXX_FLAGS_RELEASE_INIT "-MD -Ox -DNDEBUG")
  set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-MD -Os -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-MD -O2 -Zi -DNDEBUG")

  # Extra modes
  set(CMAKE_CXX_FLAGS_TESTRELEASE_INIT "-MDd -Zi -G4DEBUG_VERBOSE")
  set(CMAKE_CXX_FLAGS_MAINTAINER_INIT "-MDd -Zi")
endif()

#-----------------------------------------------------------------------
# Intel C++ Compilers - all (?) platforms
#
# Sufficient id on all platforms?
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  # Warnings
  set(CMAKE_CXX_FLAGS_INIT "-W -Wall -pedantic -Wno-non-virtual-dtor -Wno-long-long -Wwrite-strings -Wpointer-arith -Woverloaded-virtual -Wno-variadic-macros -Wshadow -fp-model precise -diag-disable=10441")
  # Use pipes rather than temp files
  set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -pipe")

  set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g")
  set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-Os -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g -DNDEBUG")

  # Extra modes
  set(CMAKE_CXX_FLAGS_TESTRELEASE_INIT "-g -G4DEBUG_VERBOSE")
  set(CMAKE_CXX_FLAGS_MAINTAINER_INIT "-g")

  # - Multithreading
  __configure_tls_models()
  set(GEANT4_MULTITHREADED_CXX_FLAGS "-pthread")

  # Linker flags
  set(CMAKE_EXE_LINKER_FLAGS "-limf")
endif()
