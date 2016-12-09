# - Geant4MakeRules_cxx
# Sets the default make rules for a CXX build, specifically the
# initialization of the compiler flags on a platform and compiler
# dependent basis
#
#
#-----------------------------------------------------------------------
# function __configure_tls_models()
#          Set available thread local storage models. Valid for GNU,
#          Clang and Intel compilers.
#
function(__configure_tls_models)
  # available models, default first
  set(_TLSMODELS initial-exec local-exec global-dynamic local-dynamic)

  set(TLSMODEL_IS_AVAILABLE ${_TLSMODELS} PARENT_SCOPE)
  foreach(_s ${_TLSMODELS})
    set(${_s}_FLAGS "-ftls-model=${_s}" PARENT_SCOPE)
  endforeach()
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

  # Additional per-mode flags
  set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g -DG4FPE_DEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -DNDEBUG")
  # Assist auto-vectorization
  set(CMAKE_CXX_FLAGS_RELEASE_INIT "${CMAKE_CXX_FLAGS_RELEASE_INIT} -fno-trapping-math -ftree-vectorize -fno-math-errno")
  set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-Os -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")

  # Remove superfluous "unused argument" "warnings" from Clang
  if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
    set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -Qunused-arguments")
  endif()

  # Though it should be the default, always use libc++ with AppleClang
  if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -stdlib=libc++")
  endif()

  # Extra Geant4 modes
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
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-MD -O2 -Zi")

  # Extra modes
  set(CMAKE_CXX_FLAGS_TESTRELEASE_INIT "-MDd -Zi -G4DEBUG_VERBOSE")
  set(CMAKE_CXX_FLAGS_MAINTAINER_INIT "-MDd -Zi")

endif()

#-----------------------------------------------------------------------
# Intel C++ Compilers - all (?) platforms
#
# Sufficient id on all platforms?
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS_INIT "-w1 -Wno-non-virtual-dtor -Wpointer-arith -Wwrite-strings -fp-model precise")
  set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g")
  set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-Os -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")

  # Extra modes
  set(CMAKE_CXX_FLAGS_TESTRELEASE_INIT "-g -G4DEBUG_VERBOSE")
  set(CMAKE_CXX_FLAGS_MAINTAINER_INIT "-g")

  # - Multithreading
  __configure_tls_models()
  set(GEANT4_MULTITHREADED_CXX_FLAGS "-pthread")

  # Linker flags
  set(CMAKE_EXE_LINKER_FLAGS "-limf")
endif()


#-----------------------------------------------------------------------
# Ye Olde *NIX/Compiler Systems
# NB: *NOT* Supported... Only provided as legacy.
# None are tested...
# Whilst these use flags taken from existing Geant4 setup, may want to see if
# CMake defaults on these platforms are good enough...
#
if(UNIX AND NOT CMAKE_COMPILER_IS_GNUCXX)
  #---------------------------------------------------------------------
  # IBM xlC compiler
  #
  if(CMAKE_CXX_COMPILER MATCHES "xlC")
    set(CMAKE_CXX_FLAGS_INIT "")
    set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g -qdbextra -qcheck=all -qfullpath -qtwolink -+")
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O3 -qtwolink -+")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-O3 -qtwolink -+")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O3 -g -qdbextra -qcheck=all -qfullpath -qtwolink -+")
  endif()

  #---------------------------------------------------------------------
  # HP aC++ Compiler
  #
  if(CMAKE_CXX_COMPILER MATCHES "aCC")
    set(CMAKE_CXX_FLAGS_INIT "+DAportable +W823")
    set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g")
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "+O2 +Onolimit")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-O2 +Onolimit")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 +Onolimit -g")
  endif()

  #---------------------------------------------------------------------
  # IRIX MIPSpro CC Compiler
  #
  if(CMAKE_CXX_COMPILER MATCHES "CC" AND CMAKE_SYSTEM_NAME MATCHES "IRIX")
    set(CMAKE_CXX_FLAGS_INIT "-ptused -DSOCKET_IRIX_SOLARIS")
    set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g")
    set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O -OPT:Olimit=5000")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-O -OPT:Olimit=5000")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O -OPT:Olimit=5000 -g")
  endif()

  #---------------------------------------------------------------------
  # SunOS CC Compiler
  # - CMake may do a reasonable job on its own here...
endif()

