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
include(IntelCompileFeatures)

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
# Provide optional file level parallelization with MSVC compiler
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
# Configure/Select C++ Standard
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

# - GEANT4_BUILD_CXXSTD
# Choose C++ Standard to build against from supported list. Allow user
# to supply it as a simple year or as 'c++XY'. If the latter, post process
# to remove the 'c++'
# Mark as advanced because most users will not need it
enum_option(GEANT4_BUILD_CXXSTD
  DOC "C++ Standard to compile against"
  VALUES 11 14 c++11 c++14
  CASE_INSENSITIVE
  )

string(REGEX REPLACE "^c\\+\\+" "" GEANT4_BUILD_CXXSTD "${GEANT4_BUILD_CXXSTD}")
mark_as_advanced(GEANT4_BUILD_CXXSTD)
geant4_add_feature(GEANT4_BUILD_CXXSTD "Compiling against C++ Standard '${GEANT4_BUILD_CXXSTD}'")

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
# config files. Needed to support late CMake 2.8 where compile features
# are not available.
set(GEANT4_CXXSTD_FLAGS "${CMAKE_CXX${GEANT4_BUILD_CXXSTD}_STANDARD_COMPILE_OPTION}")

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

# On WIN32, we need to build the genwindef application to create export
# def files for building DLLs.
# We only use it as a helper application at the moment so we exclude it from
# the ALL target.
# TODO: We could move this section into the Geant4MacroLibraryTargets.cmake
# if it can be protected so that the genwindef target wouldn't be defined
# more than once... Put it here for now...
if(WIN32)
  # Assume the sources are co-located
  get_filename_component(_genwindef_src_dir ${CMAKE_CURRENT_LIST_FILE} PATH)
  add_executable(genwindef EXCLUDE_FROM_ALL
    ${_genwindef_src_dir}/genwindef/genwindef.cpp
    ${_genwindef_src_dir}/genwindef/LibSymbolInfo.h
    ${_genwindef_src_dir}/genwindef/LibSymbolInfo.cpp)
endif()

#------------------------------------------------------------------------
# Setup symbol visibility (library interface)
# We need to define that we're building Geant4
#

#------------------------------------------------------------------------
# Optional build of examples - only intended for testing
#
option(GEANT4_BUILD_EXAMPLES "Build all the examples of the project" OFF)
GEANT4_ADD_FEATURE(GEANT4_BUILD_EXAMPLES "Build all the examples of the project")
mark_as_advanced(GEANT4_BUILD_EXAMPLES)


#-----------------------------------------------------------------------
# Integration and unit tests
# - "ENABLE_TESTING" means all tests under tests/
option(GEANT4_ENABLE_TESTING "Enable and define all the tests of the project" OFF)
GEANT4_ADD_FEATURE(GEANT4_ENABLE_TESTING "Enable and define all the tests of the project")
mark_as_advanced(GEANT4_ENABLE_TESTING)

# - "BUILD_TESTS" means all 'tests' in individual categories.
option(GEANT4_BUILD_TESTS "Build all the tests of the project" OFF)
GEANT4_ADD_FEATURE(GEANT4_BUILD_TESTS "Build all the tests of the project")
mark_as_advanced(GEANT4_BUILD_TESTS)


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
elseif(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
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
# - Check that install dirs have been defined as we want to match the
#   output layout!
if((NOT DEFINED CMAKE_INSTALL_BINDIR) OR (NOT DEFINED CMAKE_INSTALL_LIBDIR))
  message(FATAL_ERROR "Cannot configure build output dirs as install directories have not yet been defined")
endif()

# - Single root dir of all products
set(BASE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/BuildProducts")

# - Default outputs for different products, will be used by single mode
#   generators. Creates the structure:
#
# BuildProducts/
# +- bin/
# +- lib/
#    +- <GEANT4_SYSTEM>-<GEANT4_COMPILER> -> symlink -> .
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${BASE_OUTPUT_DIRECTORY}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${BASE_OUTPUT_DIRECTORY}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${BASE_OUTPUT_DIRECTORY}/${CMAKE_INSTALL_LIBDIR}")

# - Create libdir/softlink to fool geant4make, but only for single mode case
if(UNIX AND NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT EXISTS "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${GEANT4_SYSTEM}-${GEANT4_COMPILER}")
    file(MAKE_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink . ${GEANT4_SYSTEM}-${GEANT4_COMPILER}
      WORKING_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
      )
  endif()
endif()

# - For multiconfig generators, we create the same structure once for each
#   mode. Results in the structure:
#
# BuildProducts/
# +- Release/
# |  +- bin/
# |  +- lib/
# +- Debug/
# |  +- bin/
# |  +- lib/
# | ...
#
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

