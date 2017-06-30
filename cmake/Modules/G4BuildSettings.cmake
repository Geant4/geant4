#.rst:
# G4BuildSettings.cmake
# ---------------------
#
# Set defaults for core CMake settings for compiling and linking
# Geant4 libraries and tests,
#
# In addition to the core compiler/linker flags (configured in the
# Geant4MakeRules_<LANG>.cmake files) for Geant4, the build may require
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
  if(NOT CMAKE_BUILD_TYPE)
    # Default to a Release build if nothing else...
    set(CMAKE_BUILD_TYPE Release
      CACHE STRING "Choose the type of build, options are: None Release TestRelease MinSizeRel Debug RelWithDebInfo MinSizeRel Maintainer."
      FORCE
      )
  else()
    # Force to the cache, but use existing value.
    set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}"
      CACHE STRING "Choose the type of build, options are: None Release TestRelease MinSizeRel Debug RelWithDebInfo MinSizeRel Maintainer."
      FORCE
      )
  endif()
else()
  # Multimode tools like VS, Xcode
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
# - ``GEANT4_BUILD_CXXSTD`` (Allowed values: 11, 14, c++11, c++14)
#
#   - Choose C++ Standard to build against from supported list.
#
enum_option(GEANT4_BUILD_CXXSTD
  DOC "C++ Standard to compile against"
  VALUES 11 14 c++11 c++14
  CASE_INSENSITIVE
  )

string(REGEX REPLACE "^c\\+\\+" "" GEANT4_BUILD_CXXSTD "${GEANT4_BUILD_CXXSTD}")
mark_as_advanced(GEANT4_BUILD_CXXSTD)

# Require at least C++11 with no extensions and the following features
set(CMAKE_CXX_EXTENSIONS OFF)
set(GEANT4_TARGET_COMPILE_FEATURES
  cxx_alias_templates
  cxx_auto_type
  cxx_delegating_constructors
  cxx_enum_forward_declarations
  cxx_explicit_conversions
  cxx_final
  cxx_lambdas
  cxx_nullptr
  cxx_override
  cxx_range_for
  cxx_strong_enums
  cxx_uniform_initialization
  # Features that MSVC 18.0 cannot support but in list of Geant4 coding
  # guidelines - to be required once support for that compiler is dropped.
  # Version 10.2 is coded without these being required.
  #cxx_deleted_functions
  #cxx_generalized_initializers
  #cxx_constexpr
  #cxx_inheriting_constructors
  )

# If a standard higher than 11 has been selected, check that compiler has
# at least one feature from that standard and append these to the required
# feature list
if(GEANT4_BUILD_CXXSTD GREATER 11)
  if(CMAKE_CXX${GEANT4_BUILD_CXXSTD}_COMPILE_FEATURES)
    list(APPEND GEANT4_TARGET_COMPILE_FEATURES ${CMAKE_CXX${GEANT4_BUILD_CXXSTD}_COMPILE_FEATURES})
  else()
    message(FATAL_ERROR "Geant4 requested to be compiled against C++ standard '${GEANT4_BUILD_CXXSTD}'\nbut detected compiler '${CMAKE_CXX_COMPILER_ID}', version '${CMAKE_CXX_COMPILER_VERSION}'\ndoes not support any features of that standard")
  endif()
endif()

# - Check for Standard Library Implementation Features
# Smart pointers are a library implementation feature
# Hashed containers are a library implementation feature
# Random numbers are a library implementation feature?
# - Thread local? Yes, though on AppleClang platforms, see this:
  #http://stackoverflow.com/questions/28094794/why-does-apple-clang-disallow-c11-thread-local-when-official-clang-supports
# An example of where a workaround is needed
# Rest of concurrency a library implementation feature

# Add Definition to flags for temporary back compatibility
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4USE_STD11")

# Hold any appropriate compile flag(s) in variable for later export to
# config files. Needed to support clients using late CMake 2.8 where compile features
# are not available.
set(GEANT4_CXXSTD_FLAGS "${CMAKE_CXX${GEANT4_BUILD_CXXSTD}_STANDARD_COMPILE_OPTION}")

#-----------------------------------------------------------------------
# Multithreading
#-----------------------------------------------------------------------
#.rst:
# - ``GEANT4_BUILD_MULTITHREADED`` (Default: ``OFF``)
#
#   If set to ``ON``, compile Geant4 with support for multithreading.
#   At present, this is only supported on Unix platforms.
#
option(GEANT4_BUILD_MULTITHREADED "Enable multithreading in Geant4" OFF)

if(WIN32)
  mark_as_advanced(GEANT4_BUILD_MULTITHREADED)
endif()

geant4_add_feature(GEANT4_BUILD_MULTITHREADED "Build multithread enabled libraries")

#.rst:
# - ``GEANT4_BUILD_TLS_MODEL`` (Default: ``initial-exec``)
#
#   If ``GEANT4_BUILD_MULTITHREADED`` is ``ON`` and a GNU/Clang/Intel
#   compiler is being used, then this option may be used to choose the
#   Thread Local Storage model. A default of ``initial-exec`` is used
#   to provide optimize performance for applications directly linking
#   to the Geant4 libraries
#
if(GEANT4_BUILD_MULTITHREADED)
  # - Need Thread Local Storage support (POSIX)
  if(UNIX)
    check_cxx_source_compiles("__thread int i; int main(){return 0;}" HAVE_TLS)
    if(NOT HAVE_TLS)
      message(FATAL_ERROR "Configured compiler ${CMAKE_CXX_COMPILER} does not support thread local storage")
    endif()
  endif()

  # - Emit warning on Windows - message will format oddly on CMake prior
  # to 2.8, but still print
  if(WIN32)
    message(WARNING "GEANT4_BUILD_MULTITHREADED IS NOT SUPPORTED on Win32. This option should only be activated by developers")
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
    geant4_add_feature(GEANT4_BUILD_TLS_MODEL "Building with TLS model '${GEANT4_BUILD_TLS_MODEL}'")

    set(GEANT4_MULTITHREADED_CXX_FLAGS "${GEANT4_MULTITHREADED_CXX_FLAGS} ${${GEANT4_BUILD_TLS_MODEL}_FLAGS}")
  endif()

  # Set Defs/Compiler Flags
  # TODO: Migrate def to header
  add_definitions(-DG4MULTITHREADED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GEANT4_MULTITHREADED_CXX_FLAGS}")
endif()

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

# TODO: Migrate this to header
if(GEANT4_BUILD_STORE_TRAJECTORY)
  add_definitions(-DG4_STORE_TRAJECTORY)
endif()

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

# TODO: Migrate this to header
if(GEANT4_BUILD_VERBOSE_CODE)
  add_definitions(-DG4VERBOSE)
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

#.rst:
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

# On WIN32, we need to build the genwindef application to create export
# def files for building DLLs.
# We only use it as a helper application at the moment so we exclude it from
# the ALL target.
# TODO: We could move this section into the Geant4MacroLibraryTargets.cmake
# if it can be protected so that the genwindef target wouldn't be defined
# more than once... Put it here for now...
# TODO: Investigate use of CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS from CMake 3.4
#       and/or other forms of symbol visibility
if(WIN32)
  # Assume the sources are co-located
  get_filename_component(_genwindef_src_dir ${CMAKE_CURRENT_LIST_FILE} PATH)
  add_executable(genwindef EXCLUDE_FROM_ALL
    ${_genwindef_src_dir}/genwindef/genwindef.cpp
    ${_genwindef_src_dir}/genwindef/LibSymbolInfo.h
    ${_genwindef_src_dir}/genwindef/LibSymbolInfo.cpp)
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


