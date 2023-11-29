#.rst:
# G4BuildSettings.cmake
# ---------------------
#
# Set defaults for core CMake settings for compiling and linking
# Geant4 libraries and tests,
#
# In addition to the core compiler/linker flags (configured in the
# G4MakeRules_<LANG>.cmake files) for Geant4, the build may require
# further configuration. This module performs this task whicj includes:
#
#  1) Extra build modes for developers
#  2) Additional compiler definitions to assist visualization or optimize
#     performance.
#  3) Additional compiler flags which may be added optionally.
#  4) Whether to build shared and/or static libraries.
#  5) Whether to build libraries in global or granular format.
#

#-----------------------------------------------------------------
# License and Disclaimer
#
# The  Geant4 software  is  copyright of the Copyright Holders  of
# the Geant4 Collaboration.  It is provided  under  the terms  and
# conditions of the Geant4 Software License,  included in the file
# LICENSE and available at  http://cern.ch/geant4/license .  These
# include a list of copyright holders.
#
# Neither the authors of this software system, nor their employing
# institutes,nor the agencies providing financial support for this
# work  make  any representation or  warranty, express or implied,
# regarding  this  software system or assume any liability for its
# use.  Please see the license in the file  LICENSE  and URL above
# for the full disclaimer and the limitation of liability.
#
# This  code  implementation is the result of  the  scientific and
# technical work of the GEANT4 collaboration.
# By using,  copying,  modifying or  distributing the software (or
# any work based  on the software)  you  agree  to acknowledge its
# use  in  resulting  scientific  publications,  and indicate your
# acceptance of all terms of the Geant4 Software license.
#
#-----------------------------------------------------------------

if(NOT __G4BUILDSETTINGS_INCLUDED)
  set(__G4BUILDSETTINGS_INCLUDED TRUE)
else()
  return()
endif()

#-----------------------------------------------------------------------
#.rst:
# Included Modules
# ^^^^^^^^^^^^^^^^
#
# This module includes the following modules:
#
include(CheckCXXSourceCompiles)
include(IntelCompileFeatures)
include(MSVCCompileFeatures)

#-----------------------------------------------------------------------
#.rst:
# Build Modes
# ^^^^^^^^^^^
#
# Geant4 supports the standard CMake build types/configurations of
#
# - ``Release``
# - ``Debug``
# - ``MinSizeRel``
# - ``RelWithDebInfo``
#
# For single build type tools like `make` and `ninja`, the build
# defaults to ``Release`` mode unless ``CMAKE_BUILD_TYPE`` is set.
# Multi build type tools like Visual Studio and Xcode will default
# to
#
# On non-WIN32 platforms, two types specifically for development are added:
#
if(NOT WIN32)

#.rst:
# - ``Debug_FPE``
#   For debugging with full Floating Point Exception checking
#
set(CMAKE_CXX_FLAGS_DEBUG_FPE "${CMAKE_CXX_FLAGS_DEBUG_FPE_INIT}"
  CACHE STRING "Flags used by the compiler during Debug_FPE builds"
  )
mark_as_advanced(CMAKE_CXX_FLAGS_DEBUG_FPE)

#.rst:
# - ``TestRelease``:
#   For trial production and extended testing. It has verbose
#   output, has debugging symbols, and adds definitions to allow FPE
#   and physics conservation law testing where supported.
#
set(CMAKE_CXX_FLAGS_TESTRELEASE "${CMAKE_CXX_FLAGS_TESTRELEASE_INIT}"
  CACHE STRING "Flags used by the compiler during TestRelease builds"
  )
mark_as_advanced(CMAKE_CXX_FLAGS_TESTRELEASE)

#.rst:
# - ``Maintainer``:
#   For development of the toolkit. It adds debugging, and enables the use
#   of library specific debugging via standardized definitions.
#
set(CMAKE_CXX_FLAGS_MAINTAINER "${CMAKE_CXX_FLAGS_MAINTAINER_INIT}"
  CACHE STRING "Flags used by the compiler during Maintainer builds"
  )
mark_as_advanced(CMAKE_CXX_FLAGS_MAINTAINER)

#.rst:
# Compiler flags specific to these build types are set in the cache, and
# the types are added to the ``CMAKE_BUILD_TYPES`` cache string and to
# ``CMAKE_CONFIGURATION_TYPES`` if appropriate to the build tool being used.
#
if(NOT CMAKE_CONFIGURATION_TYPES)
  # Single mode build tools like Make, Ninja,
  set(__g4buildmodes "" Release TestRelease MinSizeRel Debug Debug_FPE RelWithDebInfo MinSizeRel Maintainer)
  if(NOT CMAKE_BUILD_TYPE)
    # Default to a Release build if nothing else...
    set(CMAKE_BUILD_TYPE Release
      CACHE STRING "Choose the type of build, options are: None Release TestRelease MinSizeRel Debug Debug_FPE RelWithDebInfo MinSizeRel Maintainer."
      FORCE
      )
  else()
    # Force to the cache, but use existing value.
    set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}"
      CACHE STRING "Choose the type of build, options are: None Release TestRelease MinSizeRel Debug Debug_FPE RelWithDebInfo MinSizeRel Maintainer."
      FORCE
      )
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${__g4buildmodes})
else()
  # Multimode tools like VS, Xcode
  list(APPEND CMAKE_CONFIGURATION_TYPES Debug_FPE)
  list(APPEND CMAKE_CONFIGURATION_TYPES TestRelease)
  list(APPEND CMAKE_CONFIGURATION_TYPES Maintainer)
  list(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}"
    CACHE STRING "Geant4 configurations for multimode build tools"
    FORCE
    )
endif()

endif() #NOT WIN32


#-----------------------------------------------------------------------
#.rst:
# Compilation Options
# ^^^^^^^^^^^^^^^^^^^
#
# The following options affect the compilation of the Geant4 libraries
#

#.rst
# - ``CMAKE_CXX_STANDARD`` (Allowed values: 17, 20, 23)
# - ``GEANT4_BUILD_CXXSTD`` (DEPRECATED)
#
#   - Choose C++ Standard to build against from supported list.
#   - C++23 awareness is only available from CMake 3.20 onwards
set(__g4_default_cxxstd 17 20)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.20)
  list(APPEND __g4_default_cxxstd 23)
endif()

if(GEANT4_BUILD_CXXSTD)
  message(WARNING
    "The GEANT4_BUILD_CXXSTD argument is deprecated.\n"
    "Please use CMAKE_CXX_STANDARD to set the required C++ Standard. "
    "The value supplied will be used to set CMAKE_CXX_STANDARD, but this behaviour will "
    "be removed in future releases.")
  set(CMAKE_CXX_STANDARD ${GEANT4_BUILD_CXXSTD})
endif()

enum_option(CMAKE_CXX_STANDARD
  DOC "C++ Standard to compile against"
  VALUES ${__g4_default_cxxstd}
  CASE_INSENSITIVE
  )

mark_as_advanced(CMAKE_CXX_STANDARD)
geant4_add_feature(CMAKE_CXX_STANDARD "Compiling against C++ Standard '${CMAKE_CXX_STANDARD}'")

# Setup required CXX variables: must not use vendor extensions, must require standard
# Set convenience variable for compile features to use
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(GEANT4_TARGET_COMPILE_FEATURES cxx_std_${CMAKE_CXX_STANDARD})

# Emit early fatal error if we don't have a record of the requested standard
# in CMake's known features. CMake *will* check this for us as well, but doing
# the check here only emits one error. We can also give a better hint on the
# cause/resolution
if(NOT ("${CMAKE_CXX_COMPILE_FEATURES}" MATCHES "cxx_std_${CMAKE_CXX_STANDARD}"))
  message(FATAL_ERROR
    "Geant4 requested compilation using C++ standard '${CMAKE_CXX_STANDARD}' with compiler\n"
    "'${CMAKE_CXX_COMPILER_ID}', version '${CMAKE_CXX_COMPILER_VERSION}'\n"
    "but CMake ${CMAKE_VERSION} is not aware of any support for that standard by this compiler. You may need a newer CMake and/or compiler.\n")
endif()

# - Check for Language/Standard Library Implementation Features
#
# - A very dumb wrapper round try_compile to give simple reporting
function(check_cxx_feature _flag _bindir _src)
  try_compile(${_flag} ${_bindir} ${_src} ${ARGN})
  # Can be simplified after 3.17 with CHECK_PASS/FAIL
  set(__preamble "Checking C++ feature ${_flag} -")
  if(${_flag})
    message(STATUS "${__preamble} Success")
  else()
    message(STATUS "${__preamble} Failed")
  endif()
endfunction()

# - Filesystem may be experimental/filesystem, and in stdc++fs (GCC/libstdc++), c++fs (Clang/libc++)
#   G4global will provide this as a PUBLIC dependency, if required
#   - native stdlib (nothing needed)
#   - GNU: stdc++fs
#   - Clang: c++fs
check_cxx_feature(CXXSTDLIB_FILESYSTEM_NATIVE
  ${PROJECT_BINARY_DIR}/cxx_filesystem/native
  ${PROJECT_SOURCE_DIR}/cmake/Modules/check_cxx_filesystem.cc)

if(NOT CXXSTDLIB_FILESYSTEM_NATIVE)
  # GNU libstdc++fs
  check_cxx_feature(CXXSTDLIB_FILESYSTEM_STDCXXFS
    ${PROJECT_BINARY_DIR}/cxx_filesystem/stdc++fs
    ${PROJECT_SOURCE_DIR}/cmake/Modules/check_cxx_filesystem.cc
    LINK_LIBRARIES stdc++fs)
  # LLVM libc++fs
  check_cxx_feature(CXXSTDLIB_FILESYSTEM_CXXFS
    ${PROJECT_BINARY_DIR}/cxx_filesystem/c++fs
    ${PROJECT_SOURCE_DIR}/cmake/Modules/check_cxx_filesystem.cc
    LINK_LIBRARIES c++fs)
endif()
# Derive name of library to link to
# Default is native...
if(CXXSTDLIB_FILESYSTEM_NATIVE)
  set(GEANT4_CXX_FILESYSTEM_LIBRARY )
elseif(CXXSTDLIB_FILESYSTEM_STDCXXFS)
  set(GEANT4_CXX_FILESYSTEM_LIBRARY "stdc++fs")
elseif(CXXSTDLIB_FILESYSTEM_CXXFS)
  set(GEANT4_CXX_FILESYSTEM_LIBRARY "c++fs")
else()
  message(FATAL_ERROR "No support for C++ filesystem found for compiler '${CMAKE_CXX_COMPILER_ID}', '${CMAKE_CXX_COMPILER_VERSION}'")
endif()

# Hold any appropriate compile flag(s) in variable for later export to
# non-cmake config files
set(GEANT4_CXXSTD_FLAGS "${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION}")

#.rst:
# - ``GEANT4_BUILD_SANITIZER`` (Allowed values: none address thread undefined)
#
#   - Build libraries with instrumentation of the given sanitizer
#   - Only for GNU/Clang at the moment
#   - VS 2019 16.9 supports asan: https://devblogs.microsoft.com/cppblog/address-sanitizer-for-msvc-now-generally-available/)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|.*Clang")
  enum_option(GEANT4_BUILD_SANITIZER
    DOC "Build libraries with sanitizer instrumentation"
    VALUES none address thread undefined
    CASE_INSENSITIVE)

  mark_as_advanced(GEANT4_BUILD_SANITIZER)

  if(NOT (GEANT4_BUILD_SANITIZER STREQUAL "none"))
    # Add flags - longer term, make compile/link options
    # frame pointer flag to get more meaningful stack traces
    # May need others for better/reliable output
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=${GEANT4_BUILD_SANITIZER}")

    geant4_add_feature(GEANT4_BUILD_SANITIZER "Compiling/linking with sanitizer '${GEANT4_BUILD_SANITIZER}'")
  endif()
endif()

#.rst:
# - ``GEANT4_BUILD_ENABLE_ASSERTIONS``
#   - Build libraries with assertions enabled even in release build types (Release, RelWithDebInfo, MinSizeRel)
#
option(GEANT4_BUILD_ENABLE_ASSERTIONS "Enable assertions regardless of build mode" OFF)
geant4_add_feature(GEANT4_BUILD_ENABLE_ASSERTIONS "Compiling Geant4 with assertions enabled in all build types")
mark_as_advanced(GEANT4_BUILD_ENABLE_ASSERTIONS)

# "Dumb" strip of NDEBUG definition from build type flags
if(GEANT4_BUILD_ENABLE_ASSERTIONS)
  foreach(__build_type Release MinSizeRel RelWithDebInfo)
    string(TOUPPER ${__build_type} __build_type_uc)
    foreach(__lang C CXX)
      set(FLAGSET_TO_MODIFY "CMAKE_${__lang}_FLAGS_${__build_type_uc}")
      string(REGEX REPLACE "(^| )[/-]D *NDEBUG($| )" " " new_flags "${${FLAGSET_TO_MODIFY}}")
      set(${FLAGSET_TO_MODIFY} "${new_flags}")
    endforeach()
  endforeach()
endif()

# Copy .clang-tidy file to build dir
# - clang-tidy/run-clang-tidy state that by default:
#
#   "clang-tidy will attempt to find a file named .clang-tidy for each source file in its parent directories."
#
#   Whilst this does appear to be the case, `run-clang-tidy` will report enabled checks
#   as the default set, or that of any .clang-tidy found in parent directories of the build dir
#   To avoid confusion, copy .clang-tidy into build directory so that run-clang-tidy
#   should report correctly
#
if(EXISTS "${PROJECT_SOURCE_DIR}/.clang-tidy")
  configure_file("${PROJECT_SOURCE_DIR}/.clang-tidy" "${PROJECT_BINARY_DIR}/.clang-tidy" COPYONLY)
endif()

#-----------------------------------------------------------------------
# Multithreading
#-----------------------------------------------------------------------
#.rst:
# - ``GEANT4_BUILD_MULTITHREADED`` (Default: ``OFF``)
#
#   If set to ``ON``, compile Geant4 with support for multithreading.
#   On Win32, this option requires use of static libraries only
#
option(GEANT4_BUILD_MULTITHREADED "Enable multithreading in Geant4" ON)
geant4_add_feature(GEANT4_BUILD_MULTITHREADED "Build multithread enabled libraries")

#.rst:
# - ``GEANT4_BUILD_TLS_MODEL`` (Default: ``initial-exec``)
#
#   If ``GEANT4_BUILD_MULTITHREADED`` is ``ON`` and a GNU/Clang/Intel
#   compiler is being used, then this option may be used to choose the
#   Thread Local Storage model. A default of ``initial-exec`` is used
#   to provide optimal performance for applications directly linking
#   to the Geant4 libraries. The dummy ``auto`` model may be selected
#   to leave model selection to the compiler, with no additional flag(s)
#   passed to the compiler.
#
if(GEANT4_BUILD_MULTITHREADED)
  # - Need Thread Local Storage support (POSIX)
  if(UNIX)
    check_cxx_source_compiles("__thread int i; int main(){return 0;}" HAVE_TLS)
    if(NOT HAVE_TLS)
      message(FATAL_ERROR "Configured compiler ${CMAKE_CXX_COMPILER} does not support thread local storage")
    endif()
  endif()

  # - Allow advanced users to select the thread local storage model,
  # if the compiler supports it, defaulting to that recommended by Geant4
  if(TLSMODEL_IS_AVAILABLE)
    enum_option(GEANT4_BUILD_TLS_MODEL
      DOC "Build libraries with Thread Local Storage model"
      VALUES ${TLSMODEL_IS_AVAILABLE}
      CASE_INSENSITIVE
    )
    mark_as_advanced(GEANT4_BUILD_TLS_MODEL)

    if(GEANT4_BUILD_TLS_MODEL STREQUAL "auto")
      geant4_add_feature(GEANT4_BUILD_TLS_MODEL "Building without explicit TLS model")
    else()
      geant4_add_feature(GEANT4_BUILD_TLS_MODEL "Building with TLS model '${GEANT4_BUILD_TLS_MODEL}'")
    endif()

    set(GEANT4_MULTITHREADED_CXX_FLAGS "${GEANT4_MULTITHREADED_CXX_FLAGS} ${${GEANT4_BUILD_TLS_MODEL}_TLSMODEL_FLAGS}")
  endif()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GEANT4_MULTITHREADED_CXX_FLAGS}")
endif()

#-----------------------------------------------------------------------
# Physics
#-----------------------------------------------------------------------
#.rst:
# - ``GEANT4_BUILD_PHP_AS_HP`` (Default: OFF)
#
#   - Build ParticleHP as HP.
#
set(_default_build_php_as_hp OFF)
option(GEANT4_BUILD_PHP_AS_HP "Build ParticleHP as HP" ${_default_build_php_as_hp})
mark_as_advanced(GEANT4_BUILD_PHP_AS_HP)
geant4_add_feature(GEANT4_BUILD_PHP_AS_HP "Building ParticleHP as HP")

#-----------------------------------------------------------------------
# Miscellaneous
#-----------------------------------------------------------------------
#.rst:
# - ``GEANT4_BUILD_STORE_TRAJECTORY`` (Default: ON)
#
#   - Switching off can improve performance. Needs to be on
#     for visualization to work fully. Mark as advanced because most users
#     should not need to worry about it.
#
option(GEANT4_BUILD_STORE_TRAJECTORY
  "Store trajectories in event processing. Switch off for improved performance but note that visualization of trajectories will not be possible"
  ON)
mark_as_advanced(GEANT4_BUILD_STORE_TRAJECTORY)

#.rst:
# - ``GEANT4_BUILD_VERBOSE_CODE`` (Default: ON)
#
#   - Switching off can improve performance, but at the cost
#     of fewer informational or warning messages.
#     Mark as advanced because most users should not need to worry about it.
#
option(GEANT4_BUILD_VERBOSE_CODE
  "Enable verbose output from Geant4 code. Switch off for better performance at the cost of fewer informational messages or warnings"
  ON)
mark_as_advanced(GEANT4_BUILD_VERBOSE_CODE)

#.rst:
# - ``GEANT4_BUILD_BUILTIN_BACKTRACE`` (Unix only, Default: OFF)
#
#   - Setting to ``ON`` will build in automatic signal handling for ``G4RunManager``
#     through ``G4Backtrace``. Applications requiring/implenting their own signal
#     handling should not enable this option.
#
if(UNIX)
  option(GEANT4_BUILD_BUILTIN_BACKTRACE
    "Enable automatic G4Backtrace signal handling in G4RunManager. Switch off for applications implementing their own signal handling"
    OFF)
  mark_as_advanced(GEANT4_BUILD_BUILTIN_BACKTRACE)
endif()

#.rst:
# - ``GEANT4_BUILD_MSVC_MP`` (Windows only, Default: OFF)
#
#    - Provide optional file level parallelization with MSVC compiler.
#

# NB: This will only work if the build tool performs compilations
# with multiple sources, e.g. "cl.exe a.cc b.cc c.cc"
if(MSVC)
  option(GEANT4_BUILD_MSVC_MP "Use /MP option with MSVC for file level parallel builds" OFF)
  mark_as_advanced(GEANT4_BUILD_MSVC_MP)

  if(GEANT4_BUILD_MSVC_MP)
    set(CMAKE_CXX_FLAGS "/MP ${CMAKE_CXX_FLAGS}")
  endif()
endif()


#-----------------------------------------------------------------------
#.rst:
# Library Mode and Link Options
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# - ``BUILD_SHARED_LIBS`` (Default: ``ON``)
#
#   - Build shared libraries (`.so` on Unix, `.dylib` on macOS, `.dll/.lib` on Windows)
#
# - ``BUILD_STATIC_LIBS`` (Default: ``OFF``)
#
#   - Build static libraries (`.a` on Unices, `.lib` on Windows)
#
# Both options are adavanced and may be selected together to build
# both types. If neither is selected, an error is emitted.
#
option(BUILD_SHARED_LIBS "Build Geant4 shared libraries" ON)
option(BUILD_STATIC_LIBS "Build Geant4 static libraries" OFF)
mark_as_advanced(BUILD_SHARED_LIBS BUILD_STATIC_LIBS)

# Because both could be switched off accidently, FATAL_ERROR if neither
# option has been selected.
if(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
  message(FATAL_ERROR "Neither static nor shared libraries will be built")
endif()

# Always build global libraries - always FATAL_ERROR if old
# granular library switch is set, e.g. from command line
if(GEANT4_BUILD_GRANULAR_LIBS)
  message(FATAL_ERROR " Granular libraries are no longer supported!")
endif()

#------------------------------------------------------------------------
#.rst:
# Build Output Directories
# ------------------------
#
# Because of the highly nested structure of Geant4, built libraries and
# programs will be distributed throughout the build tree. To simplify
# use and testing of the build, CMake is setup to output these products
# to a directory hierarchy based on the configured install locations.

# - Match output layout, so must have variables from GNUInstallDirs...
if((NOT DEFINED CMAKE_INSTALL_BINDIR) OR (NOT DEFINED CMAKE_INSTALL_LIBDIR))
  message(FATAL_ERROR "Cannot configure build output dirs as install directories have not yet been defined")
endif()

#.rst
# In all cases, build products are stored in a directory tree rooted
# in a directory named ``BuildProducts`` under the ``PROJECT_BINARY_DIRECTORY``.
#
set(BASE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/BuildProducts")

#.rst:
# For single mode build generators (make, ninja), the following
#  hierarchy is used:
#
#  .. code-block:: console
#
#    +- <PROJECT_BINARY_DIR>/
#       +- BuildProducts/
#          +- <CMAKE_INSTALL_BINDIR>/
#             +- ... "runtime" targets ...
#          +- <CMAKE_INSTALL_LIBDIR>/
#             +- ... "library" and "archive" targets ...
#

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${BASE_OUTPUT_DIRECTORY}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${BASE_OUTPUT_DIRECTORY}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${BASE_OUTPUT_DIRECTORY}/${CMAKE_INSTALL_LIBDIR}")

#.rst:
#  For multimode build generators (Xcode, Visual Studio), each mode
#  is separated using the hierarchy
#
#  .. code-block:: console
#
#    +- <PROJECT_BINARY_DIR>
#       +- BuildProducts/
#          +- <CONFIG>/
#             +- <CMAKE_INSTALL_BINDIR>/
#                +- ... "runtime" targets ...
#             +- <CMAKE_INSTALL_LIBDIR>/
#                +- ... "library" and "archive" targets ...
#          +- ...
#
#  where ``<CONFIG>`` is repeated for each build configuration listed in
#  :cmake:variable:`CMAKE_CONFIGURATION_TYPES <cmake:variable:CMAKE_CONFIGURATION_TYPES>`, e.g. Release, Debug, RelWithDebInfo etc.

foreach(_conftype ${CMAKE_CONFIGURATION_TYPES})
  string(TOUPPER ${_conftype} _conftype_uppercase)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${_conftype_uppercase}
    "${BASE_OUTPUT_DIRECTORY}/${_conftype}/${CMAKE_INSTALL_BINDIR}"
    )
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${_conftype_uppercase}
    "${BASE_OUTPUT_DIRECTORY}/${_conftype}/${CMAKE_INSTALL_LIBDIR}"
    )
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${_conftype_uppercase}
    "${BASE_OUTPUT_DIRECTORY}/${_conftype}/${CMAKE_INSTALL_LIBDIR}"
    )
endforeach()
