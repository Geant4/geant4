#------------------------------------------------------------------------------
# Top Level CMakeLists.txt for Geant4 Build
#
# 21st September 2010 Ben Morgan
#
# $Id: CMakeLists.txt,v 1.6 2010-12-08 16:36:34 bmorgan Exp $
# GEANT4 Tag $Name: not supported by cvs2svn $
#

#------------------------------------------------------------------------------
# Enforce out-of-source builds before anything else
#
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "Geant4 requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "CMakeCache.txt")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
endif()

#------------------------------------------------------------------------------
# Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 2.6)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
   ${CMAKE_SOURCE_DIR}/cmake/Modules/Geant4MakeRules_cxx.cmake)

set(CMAKE_USER_MAKE_RULES_OVERRIDE_Fortran
   ${CMAKE_SOURCE_DIR}/cmake/Modules/Geant4MakeRules_fortran.cmake)

#------------------------------------------------------------------------------
# Define the project
#
project(geant4)

#------------------------------------------------------------------------------
# Define the version of the project - can be used to generate sources,
# but just used for CPack the the moment
#
set(${PROJECT_NAME}_VERSION "9.4.4")
set(${PROJECT_NAME}_VERSION_MAJOR "9")
set(${PROJECT_NAME}_VERSION_MINOR "4")
set(${PROJECT_NAME}_VERSION_PATCH "4")

#------------------------------------------------------------------------------
# Set up path to internal set of CMake modules.
#
set(CMAKE_MODULE_PATH
    ${geant4_SOURCE_DIR}/cmake
    ${geant4_SOURCE_DIR}/cmake/Modules
    ${CMAKE_MODULE_PATH})


#------------------------------------------------------------------------------
# Set up Build types or configurations
# If further tuning of compiler flags is needed then it should be done here.
# (It can't be done in the make rules override section).
# However, exercise care when doing this not to override existing flags!!
# We don't do this on WIN32 platforms because of some teething issues
# with compiler specifics and linker flags
if(NOT WIN32)
    include(Geant4BuildModes)
endif(NOT WIN32)

# Optional define flags which are applicable globally

# G4_STORE_TRAJECTORY
# ON by default, switching off can improve performance. Needs to be on
# for visualization to work fully. Mark as advanced because most users
# should not need to worry about it.
# FIXES : Bug #1208
option(GEANT4_STORE_TRAJECTORY "Store trajectories in event processing. Switch off for improved performance but note that visualization of trajectories will not be possible" ON)
mark_as_advanced(GEANT4_STORE_TRAJECTORY)

if(GEANT4_STORE_TRAJECTORY)
    add_definitions(-DG4_STORE_TRAJECTORY)
endif()


# G4VERBOSE
# ON by default, switching off can improve performance, but at the cost of
# fewer informational or warning messages. Mark as advanced because most users
# should not need to worry about it.
option(GEANT4_VERBOSE_CODE 
    "Enable verbose output from Geant4 code. Switch off for better performance at the cost of fewer informational messages or warnings"
    ON)
mark_as_advanced(GEANT4_VERBOSE_CODE)

if(GEANT4_VERBOSE_CODE)
    add_definitions(-DG4VERBOSE)
endif()


#------------------------------------------------------------------------------
# Includes for common utilities
#
include(Geant4CPackBase)
include(FeatureSummary)
include(Geant4MacroUtilities)
include(CMakeDependentOption)


#------------------------------------------------------------------------------
# Geant4 library build options. 
# Only global or granular libraries can be built. Global libraries are the
# default, and provide an advanced option to switch to granular build.
# Only developers really need the granular build.
option(GEANT4_BUILD_GRANULAR_LIBS "Build Geant4 with granular libraries" OFF)
mark_as_advanced(GEANT4_BUILD_GRANULAR_LIBS)
GEANT4_ADD_FEATURE(GEANT4_BUILD_GRANULAR_LIBS "Build granular Geant4 libraries")

# Default to static libraries on win32 until declspec issues resolved
if(WIN32)
    option(BUILD_STATIC_LIBS "Build Geant4 static libraries" ON)
    option(BUILD_SHARED_LIBS "Build Geant4 dynamic libraries" OFF)
else(WIN32)
    option(BUILD_STATIC_LIBS "Build Geant4 static libraries" OFF)
    option(BUILD_SHARED_LIBS "Build Geant4 dynamic libraries" ON)
endif(WIN32)

#------------------------------------------------------------------------------
# Build output
# Output all binaries to single directory based on their type
#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/runtime)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/library)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/outputs/archive)


#------------------------------------------------------------------------------
# Geant4 defaults for data libraries
# Not required for build, so how best to handle?
#

#------------------------------------------------------------------------------
# Find required CLHEP package
# Use REQUIRED - we MUST have it. An EXACT version number is not used at 
# present, but can be added later as needed.
#
find_package(CLHEP 2.1.0.1 REQUIRED)


#------------------------------------------------------------------------------
# Geant4 User Interface Options
# Written in separate module for clarity
#
include(Geant4InterfaceOptions)

#------------------------------------------------------------------------------
# Geant4 Extra Modules Options
#
#------------------------------------------------------------------------------
# Optional support for GDML - requires XercesC package
#
option(GEANT4_USE_GDML "Build Geant4 with GDML support" OFF)
if(GEANT4_USE_GDML)
    find_package(XercesC REQUIRED)
endif(GEANT4_USE_GDML)

GEANT4_ADD_FEATURE(GEANT4_USE_GDML "Build Geant4 with GDML support")


#------------------------------------------------------------------------------
# Build G3ToG4 converter program
# The G3toG4 *library* will always be built, but the application requires
# a Fortran compiler AND CERNLIB, so is optional.
# Only on *NIX because converter requires CERNLIB, and Windows support for
# this is an unknown quantity at present (can always change later).
#
if(UNIX)
    option(GEANT4_USE_GEANT3TOGEANT4 "Build the Geant3-to-Geant4 converter application" OFF)
    if(GEANT4_USE_GEANT3TOGEANT4)
        enable_language(Fortran)
        #find_package(CERNLIB REQUIRED)
    endif(GEANT4_USE_GEANT3TOGEANT4)

    GEANT4_ADD_FEATURE(GEANT4_USE_GEANT3TOGEANT4 "Build the Geant3-to-Geant4 converter application")
endif()

#------------------------------------------------------------------------------
# Optional use of system zlib library
#
#option(GEANT4_USE_SYSTEM_ZLIB "Use the system's zlib library" OFF)
if(GEANT4_USE_SYSTEM_ZLIB)
    # This needs more work - use ITK's way of doing it as an example.
    find_package(ZLIB)
endif(GEANT4_USE_SYSTEM_ZLIB)

GEANT4_ADD_FEATURE(GEANT4_USE_SYSTEM_ZLIB "Use system zlib library")


#------------------------------------------------------------------------------
# Optional build of Python wrappers (Geant4Py)
# Suggested that this be brought under the wing of the CMake system.
# -> Requires Global/Shared library build (all platforms?)
#
#CMAKE_DEPENDENT_OPTION(GEANT4_WRAP_PYTHON "Wrap Geant4 classes into the Python language" OFF "GEANT4_BUILD_GLOBAL_LIBS;BUILD_SHARED_LIBS" OFF)
#if(GEANT4_WRAP_PYTHON)
    # We need Boost-Python and Python
    # NB: Watch for FindPythonLibs picking up the static version...
    # Known issue with CMake 2.6 module - see Bug #8389 and 2257 for details
    # of fix
    #    find_package(Boost REQUIRED COMPONENTS python)
    #find_package(PythonLibs REQUIRED)
    #endif(GEANT4_WRAP_PYTHON)

    #GEANT4_ADD_FEATURE(GEANT4_WRAP_PYTHON "Wrap Geant4 classes into the Python language")


#------------------------------------------------------------------------------
# Add the source subdirectory - all the source code for the libraries is under
# this.
#
add_subdirectory(source)


#------------------------------------------------------------------------------
# Generate any Use/Config files here once everything else has been processed
# e.g. "UseGeant4.cmake", "Geant4Config.cmake" - library dependencies etc
#
set(GEANT4_PREFIX ${CMAKE_INSTALL_PREFIX})
set(GEANT4_EXEC_PREFIX ${GEANT4_PREFIX})
set(GEANT4_BINDIR ${GEANT4_EXEC_PREFIX}/bin)
set(GEANT4_LIBDIR ${GEANT4_PREFIX}/lib)
set(GEANT4_DATAROOTDIR ${GEANT4_PREFIX}/share)
set(GEANT4_DATADIR ${GEANT4_DATAROOTDIR})
set(GEANT4_INCLUDEDIR ${GEANT4_PREFIX}/include)
set(GEANT4_DOCDIR ${GEANT4_DATAROOTDIR}/doc/geant4-${geant4_VERSION})

include(Geant4ToolchainBackwardCompatibility)
include(Geant4ConfigureConfigScript)


#------------------------------------------------------------------------------
# Install any extra files needed such as documentation and legacy Makefiles
#
# User example code - optional
option(GEANT4_INSTALL_EXAMPLES "Install all Geant4 examples" OFF)

if(GEANT4_INSTALL_EXAMPLES)
    install(DIRECTORY examples
        DESTINATION ${GEANT4_DATAROOTDIR}/geant4-${geant4_VERSION}
        COMPONENT Examples
        PATTERN "CVS" EXCLUDE
        PATTERN ".svn" EXCLUDE)
endif()

#
# Final output - show what's been enabled so that user knows what's
# happening - also useful for later problem solving!
#
GEANT4_PRINT_ENABLED_FEATURES()
