# - Setup of general build options for Geant4 Libraries
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
#

#-----------------------------------------------------------------------
# Load needed modules
#
include(CheckCXXSourceCompiles)

#-----------------------------------------------------------------------
# Set up Build types or configurations
# If further tuning of compiler flags is needed then it should be done here.
# (It can't be done in the make rules override section).
# However, exercise care when doing this not to override existing flags!!
# We don't do this on WIN32 platforms yet because of some teething issues
# with compiler specifics and linker flags
if(NOT WIN32)
  include(Geant4BuildModes)
endif()

#-----------------------------------------------------------------------
# Optional compiler definitions which are applicable globally
#
# - G4MULTITHREADED
# OFF by default. Switching on will enable multithreading, adding the
# G4MULTITHREADED definition globally and appending the relevant
# compiler flags to CMAKE_CXX_FLAGS
# Enabling the option allows advanced users to further select the
# thread local storage model if GNU/Clang/Intel compiler is used.
option(GEANT4_BUILD_MULTITHREADED "Enable multithreading in Geant4" OFF)

if(WIN32)
  mark_as_advanced(GEANT4_BUILD_MULTITHREADED)
endif()

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
  add_definitions(-DG4MULTITHREADED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GEANT4_MULTITHREADED_CXX_FLAGS}")
endif()

geant4_add_feature(GEANT4_BUILD_MULTITHREADED "Build multithread enabled libraries")

# - G4_STORE_TRAJECTORY
# ON by default, switching off can improve performance. Needs to be on
# for visualization to work fully. Mark as advanced because most users
# should not need to worry about it.
# FIXES : Bug #1208
option(GEANT4_BUILD_STORE_TRAJECTORY
  "Store trajectories in event processing. Switch off for improved performance but note that visualization of trajectories will not be possible"
  ON)
mark_as_advanced(GEANT4_BUILD_STORE_TRAJECTORY)

if(GEANT4_BUILD_STORE_TRAJECTORY)
  add_definitions(-DG4_STORE_TRAJECTORY)
endif()

# - G4VERBOSE
# ON by default, switching off can improve performance, but at the cost
# of fewer informational or warning messages. Mark as advanced because
# most users should not need to worry about it.
option(GEANT4_BUILD_VERBOSE_CODE
  "Enable verbose output from Geant4 code. Switch off for better performance at the cost of fewer informational messages or warnings"
  ON)
mark_as_advanced(GEANT4_BUILD_VERBOSE_CODE)

if(GEANT4_BUILD_VERBOSE_CODE)
  add_definitions(-DG4VERBOSE)
endif()

#-----------------------------------------------------------------------
# Optional compiler flags
#
# - GEANT4_BUILD_CXXSTD
# Choose C++ Standard to build against, if supported.
# Mark as advanced because most users will not need it.
if(CXXSTD_IS_AVAILABLE)
  enum_option(GEANT4_BUILD_CXXSTD
    DOC "C++ Standard to compile against"
    VALUES ${CXXSTD_IS_AVAILABLE}
    CASE_INSENSITIVE
    )
  mark_as_advanced(GEANT4_BUILD_CXXSTD)
  geant4_add_feature(GEANT4_BUILD_CXXSTD "Compiling against C++ Standard '${GEANT4_BUILD_CXXSTD}'")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${${GEANT4_BUILD_CXXSTD}_FLAGS}")

  # Get needed defs for selected standard
  # NB this needs refactoring because, for example, newer MSVC versions
  # will just support c++11, and the selection above isn't needed, nor
  # available, so we cannot check it in the if statement below.
  # Add the def(s) to flags, because it may be *required*
  if(GEANT4_BUILD_CXXSTD MATCHES "c\\+\\+11")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4USE_STD11")
  endif()
  if(GEANT4_BUILD_CXXSTD MATCHES "c\\+\\+14")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4USE_STD11")
  endif()
  if(GEANT4_BUILD_CXXSTD MATCHES "c\\+\\+1y")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4USE_STD11")
  endif()
  if(GEANT4_BUILD_CXXSTD MATCHES "c\\+\\+1z")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4USE_STD11")
  endif()
  if(GEANT4_BUILD_CXXSTD MATCHES "c\\+\\+17")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG4USE_STD11")
  endif()
endif()

#-----------------------------------------------------------------------
# Setup Library Format Option.
# Always build global libraries - always FATAL_ERROR if old
# granular library switch is set, e.g. from command line
if(GEANT4_BUILD_GRANULAR_LIBS)
  message(FATAL_ERROR " Granular libraries are no longer supported!")
endif()

#-----------------------------------------------------------------------
# Setup Shared and/or Static Library builds
# We name these options without a 'GEANT4_' prefix because they are
# really higher level CMake options.
# Default to building shared libraries, mark options as advanced because
# most user should not have to touch them.
option(BUILD_SHARED_LIBS "Build Geant4 shared libraries" ON)
option(BUILD_STATIC_LIBS "Build Geant4 static libraries" OFF)
mark_as_advanced(BUILD_SHARED_LIBS BUILD_STATIC_LIBS)

# Because both could be switched off accidently, FATAL_ERROR if neither
# option has been selected.
if(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
  message(FATAL_ERROR "Neither static nor shared libraries will be built")
endif()

# On MSVC, warn if both shared and static are built - this has duplicate
# symbol issues on VS2010 Express.
# TODO: Resolve and understand this issue...
if(MSVC)
  if(BUILD_SHARED_LIBS AND BUILD_STATIC_LIBS)
    message(WARNING " Building shared AND static libraries on VS2010 may result in link errors.
 You are welcome to try building both, but please be aware of this warning.
 Problems can be reported to the Geant4 Bugzilla system:

 http://bugzilla-geant4.kek.jp
    ")
  endif()
endif()

#------------------------------------------------------------------------
# Setup symbol visibility (library interface)
# We need to define that we're building Geant4
#

#------------------------------------------------------------------------
# Build options for building examples and tests
option(GEANT4_BUILD_EXAMPLES "Build all the examples of the project" OFF)
GEANT4_ADD_FEATURE(GEANT4_BUILD_EXAMPLES "Build all the examples of the project")
mark_as_advanced(GEANT4_BUILD_EXAMPLES)

option(GEANT4_BUILD_TESTS "Build all the tests of the project" OFF)
GEANT4_ADD_FEATURE(GEANT4_BUILD_TESTS "Build all the tests of the project")
mark_as_advanced(GEANT4_BUILD_TESTS)

# - Testing system only functional on 2.8 and above
if(NOT ${CMAKE_VERSION} VERSION_GREATER 2.7)
  if(GEANT4_ENABLE_TESTING)
    message(STATUS "WARNING: GEANT4_ENABLE_TESTING requires CMake >= 2.8 -- deactivating")
  endif()
  set(GEANT4_ENABLE_TESTING OFF CACHE INTERNAL "Testing NOT supported on CMake <2.8"
    FORCE)
else()
  option(GEANT4_ENABLE_TESTING "Enable and define all the tests of the project" OFF)
  GEANT4_ADD_FEATURE(GEANT4_ENABLE_TESTING "Enable and define all the tests of the project")
  mark_as_advanced(GEANT4_ENABLE_TESTING)
endif()


# On WIN32, we need to build the genwindef application to create export
# def files for building DLLs.
# We only use it as a helper application at the moment so we exclude it from
# the ALL target.
# TODO: We could move this section into the Geant4MacroLibraryTargets.cmake
# if it can be protected so that the genwindef target wouldn't be defined
# more than once... Put it here for now...
if(WIN32)
  add_definitions(-DG4LIB_BUILD_DLL)
  # Assume the sources are co-located
  get_filename_component(_genwindef_src_dir ${CMAKE_CURRENT_LIST_FILE} PATH)
  add_executable(genwindef EXCLUDE_FROM_ALL
    ${_genwindef_src_dir}/genwindef/genwindef.cpp
    ${_genwindef_src_dir}/genwindef/LibSymbolInfo.h
    ${_genwindef_src_dir}/genwindef/LibSymbolInfo.cpp)
endif()


#------------------------------------------------------------------------
# Setup Locations for Build Outputs
# Because of the highly nested structure of Geant4, targets will be
# distributed throughout this tree, potentially making usage and debugging
# difficult (especially if developers use non-CMake tools).
#
# We therefore set the output directory of runtime, library and archive
# targets to some low level directories under the build tree.
#
# On Unices, we try to make the output directory backward compatible
# with the old style 'SYSTEM-COMPILER' format so that applications may be
# built against the targets in the build tree.
#
# Note that for multi-configuration generators like VS and Xcode, these
# directories will have the configuration type (e.g. Debug) appended to
# them, so are not backward compatible with the old Make toolchain in
# these cases.
#
# Also, we only do this on UNIX because we want to avoid mixing static and
# dynamic libraries on windows until the differences are better understood.
#------------------------------------------------------------------------
# Determine the backward compatible system name
#
if(NOT WIN32)
  set(GEANT4_SYSTEM ${CMAKE_SYSTEM_NAME})
else()
  set(GEANT4_SYSTEM "WIN32")
endif()

#------------------------------------------------------------------------
# Determine the backward compatible compiler name
# NB: At present Clang detection only works on CMake > 2.8.1
if(CMAKE_COMPILER_IS_GNUCXX)
  set(GEANT4_COMPILER "g++")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(GEANT4_COMPILER "clang")

  # - Newer g++ on OS X may identify as Clang
  if(APPLE AND (CMAKE_CXX_COMPILER MATCHES ".*g\\+\\+"))
    set(GEANT4_COMPILER "g++")
  endif()

elseif(MSVC)
  set(GEANT4_COMPILER "VC")
elseif(CMAKE_CXX_COMPILER MATCHES "icpc.*|icc.*")
  set(GEANT4_COMPILER "icc")
else()
  set(GEANT4_COMPILER "UNSUPPORTED")
endif()

#-----------------------------------------------------------------------
# Set the output paths to be backward compatible on UNIX
if(NOT UNIX)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/runtime)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/library)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/archive)
else()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/runtime)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/outputs/library/${GEANT4_SYSTEM}-${GEANT4_COMPILER})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif()


