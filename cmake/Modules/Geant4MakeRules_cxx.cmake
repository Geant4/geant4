# - Geant4MakeRules_cxx
# Sets the default make rules for a CXX build, specifically the 
# initialization of the compiler flags on a platform and compiler
# dependent basis
#
#
#-----------------------------------------------------------------------
# C++ Standard Setup
#-----------------------------------------------------------------------
# If the compiler supports setting a particular C++ Standard to use
# when compiling, configure the following global variables that can
# be used later on to select the standard
#
#  CXXSTD_IS_AVAILABLE    List of C++ Standards supported by the
#                         compiler. It is empty if the compiler
#                         does not support setting the standard,
#                         e.g. MSVC.
#
#  <CXXSTD>_FLAGS         List of flags needed to compile with
#                         CXXSTD under the detected compiler.
#
# We don't add these flags to CMAKE_CXX_FLAGS_INIT as the choice
# of standard is a configure time choice, and may be changed in
# the cache. 
# 
# Settings for each compiler are handled in a dedicated function.
# Whilst we only handle GNU, Clang and Intel, these are sufficiently
# different in the required flags that individual handling is needed.
#

#-----------------------------------------------------------------------
# function __configure_cxxstd_gnu()
#          Determine version of GNU compiler and set available C++
#          Standards and flags as appropriate in the function's
#          parent scope.
#     
#          Note that this function is safe in CMake < 2.8.2 where
#          the Clang compiler is identified as the GNU compiler.
#          Clang allows the -dumpversion argument and provides
#          a value such that we will use the sensible default
#          standard of 'c++98'
#
function(__configure_cxxstd_gnu)
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
    OUTPUT_VARIABLE _gnucxx_version
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  if(_gnucxx_version VERSION_LESS 4.3)
    set(_CXXSTDS "c++98")
  elseif(_gnucxx_version VERSION_LESS 4.7)
    set(_CXXSTDS "c++98" "c++0x")
  else()
    set(_CXXSTDS "c++98" "c++0x" "c++11")
  endif()

  set(CXXSTD_IS_AVAILABLE ${_CXXSTDS} PARENT_SCOPE)
  foreach(_s ${_CXXSTDS})
    set(${_s}_FLAGS "-std=${_s}" PARENT_SCOPE)
  endforeach()
endfunction()

#-----------------------------------------------------------------------
# function __configure_cxxstd_clang()
#          Determine version of Clang compiler and set available C++
#          Standards and flags as appropriate in the function's
#          parent scope.
#
function(__configure_cxxstd_clang)
  # Hmm, Clang seems to dump -v to stderr...
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -v
    ERROR_VARIABLE _clangcxx_dumpedversion
    ERROR_STRIP_TRAILING_WHITESPACE
    )
  string(REGEX REPLACE ".*clang version ([0-9]\\.[0-9]+).*" "\\1" _clangcxx_version ${_clangcxx_dumpedversion})
  message(STATUS "Clang version : ${_clangcxx_version}")

  if(_clangcxx_version VERSION_LESS 2.9)
    set(_CXXSTDS "c++98")
  else()
    set(_CXXSTDS "c++98" "c++0x" "c++11")
  endif()

  set(CXXSTD_IS_AVAILABLE ${_CXXSTDS} PARENT_SCOPE)
  foreach(_s ${_CXXSTDS})
    set(${_s}_FLAGS "-std=${_s}" PARENT_SCOPE)
  endforeach()
endfunction()

#-----------------------------------------------------------------------
# function __configure_cxxstd_intel()
#          Determine version of Intel compiler and set available C++
#          Standards and flags as appropriate in the function's
#          parent scope.
#
function(__configure_cxxstd_intel)
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
    OUTPUT_VARIABLE _icpc_dumpedversion
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  if(_icpc_dumpedversion VERSION_LESS 11.0)
    set(_CXXSTDS "c++98")
  else()
    set(_CXXSTDS "c++98" "c++0x")
  endif()

  set(CXXSTD_IS_AVAILABLE ${_CXXSTDS} PARENT_SCOPE)
  foreach(_s ${_CXXSTDS})
    # - Intel does not support '-std=c++98'
    if(${_s} MATCHES "c\\+\\+98")
      set(${_s}_FLAGS "-ansi" PARENT_SCOPE)
    else()
      set(${_s}_FLAGS "-std=${_s}" PARENT_SCOPE)
    endif()
  endforeach()
endfunction()


#-----------------------------------------------------------------------
# DEFAULT FLAG SETTING
#-----------------------------------------------------------------------
# GNU C++ or LLVM/Clang Compiler on all(?) platforms
# NB: At present, only identifies clang correctly on CMake > 2.8.1
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS_INIT "-W -Wall -pedantic -Wno-non-virtual-dtor -Wno-long-long -Wwrite-strings -Wpointer-arith -Woverloaded-virtual -pipe")
  set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g") 
  set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-Os -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")

  # Remove superfluous "unused argument" "warnings" from Clang
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS_INIT "${CMAKE_CXX_FLAGS_INIT} -Qunused-arguments")
  endif()

  # Extra Geant4 modes
  # - TestRelease
  set(CMAKE_CXX_FLAGS_TESTRELEASE_INIT "-g -DG4DEBUG_VERBOSE -DG4FPE_DEBUG")

  # - Maintainer
  set(CMAKE_CXX_FLAGS_MAINTAINER_INIT "-g")

  # - C++ Standard Settings
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    __configure_cxxstd_clang()
  else()
    __configure_cxxstd_gnu()
  endif()
endif()


#-----------------------------------------------------------------------
# MSVC - all (?) versions
# I don't understand VS flags at all.... Guess try the CMake defaults first
# and see what happens!
if(MSVC)
  # Hmm, WIN32-VC.gmk uses dashes, but cmake uses slashes, latter probably
  # best for native build.
  set(CMAKE_CXX_FLAGS_INIT "-GR -EHsc -Zm200 -nologo -D_CONSOLE -D_WIN32 -DWIN32 -DOS -DXPNET -D_CRT_SECURE_NO_DEPRECATE")
  set(CMAKE_CXX_FLAGS_DEBUG_INIT "-MDd -Od -Zi") 
  set(CMAKE_CXX_FLAGS_RELEASE_INIT "-MD -O2 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-MD -Os -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-MD -O2 -Zi")

  # Extra modes
  set(CMAKE_CXX_FLAGS_TESTRELEASE_INIT "-MDd -Zi -G4DEBUG_VERBOSE")
  set(CMAKE_CXX_FLAGS_MAINTAINER_INIT "-MDd -Zi")

  # We may also have to set linker flags....
endif()


#-----------------------------------------------------------------------
# Intel C++ Compilers - all (?) platforms
#
# Sufficient id on all platforms?
if(CMAKE_CXX_COMPILER MATCHES "icpc.*|icc.*")
  set(CMAKE_CXX_FLAGS_INIT "-w1 -Wno-non-virtual-dtor -Wpointer-arith -Wwrite-strings -fp-model precise")
  set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g") 
  set(CMAKE_CXX_FLAGS_RELEASE_INIT "-O2 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-Os -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -g")

  # Extra modes
  set(CMAKE_CXX_FLAGS_TESTRELEASE_INIT "-g -G4DEBUG_VERBOSE")
  set(CMAKE_CXX_FLAGS_MAINTAINER_INIT "-g")

  # C++ Standard Settings
  __configure_cxxstd_intel()

  # Linker flags 
  set(CMAKE_EXE_LINKER_FLAGS "-i-dynamic -limf")
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

