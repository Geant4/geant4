#------------------------------------------------------------------------------
# Top Level CMakeLists.txt for Geant4 Build
#
# 21st September 2010 Ben Morgan
#
# $Id: CMakeLists.txt 95323 2016-02-04 15:12:37Z gcosmo $
#

#------------------------------------------------------------------------------
# - Enforce an out-of-source builds before anything else
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
# - Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

# If Policy CMP0022 exists, use OLD because Geant4 still needs
# to support versions < 2.8.12.
# This relates to LINK_INTERFACE_LIBRARIES and can be removed once
# minimum required CMake version reaches 2.8.12
if(POLICY CMP0022)
  cmake_policy(SET CMP0022 OLD)
endif()

# If Policy CMP0042 exists, use OLD to prefer the use of install names
# instead of the new @rpath default. This is temporary to suppress warnings
# on newer CMake versions.
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 OLD)
endif()

# If Policy CMP0043 exists, use OLD because Geant4 still needs
# to support Geant4 =< 2.8.12.
# This relates to the use of generator expressions to set COMPILE_DEFINITIONS
# and can be removed once the minimum required CMake version reaches 2.8.10.
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
   ${CMAKE_SOURCE_DIR}/cmake/Modules/Geant4MakeRules_cxx.cmake)

#------------------------------------------------------------------------------
# - Project definition
#
project(Geant4)


#----------------------------------------------------------------------------
# - Prepend our own CMake Modules to the search path
# This allows us to customize, add, and factor out functionality
# NB: if our custom modules include others that we don't supply, those in the
# base path will be used, so watch for incompatibilities!!
#
set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake/Modules
    ${CMAKE_MODULE_PATH})


#----------------------------------------------------------------------------
# - Add functionality provided by standard and custom modules
# See the documentation in each of these modules for further details.
#
# - Versioning. We do this here for now
set(${PROJECT_NAME}_VERSION "10.1.3")
set(${PROJECT_NAME}_VERSION_MAJOR "10")
set(${PROJECT_NAME}_VERSION_MINOR "1")
set(${PROJECT_NAME}_VERSION_PATCH "3")

# - Provide dependent options as these are needed for some Geant4 features
include(CMakeDependentOption)

# - Provide general CMake utilities for Geant4.
include(Geant4MacroUtilities)

# - Provide an 'uninstall' target.
include(CMakeUninstallTarget)

# - Provide the 'validate_sources' target.
include(Geant4ValidateSources)

# - Provide standard install directories permitting customization.
include(Geant4InstallDirs)

# - Provide options to control how Geant4 libraries are built
include(Geant4LibraryBuildOptions)

# - Provide interface to control use of optional components
include(Geant4OptionalComponents)

# - Provide interface to control use of UI/Vis components
#   Written in a separate module from other optional components because
#   there are many complex options to handle.
include(Geant4InterfaceOptions)

# - Provide options to enable wrapping of Geant4 by other languages
include(Geant4Wrapping)


#------------------------------------------------------------------------------
# Add the source and environments subdirectories
# source       : Process all the Geant4 core targets
# environments : Process optional wrappings of Geant4 (NOTYETIMPLEMENTED)
add_subdirectory(source)
#add_subdirectory(environments)


#----------------------------------------------------------------------------
# - Perform all post build tasks
# At the CMake level, this simply means that we must know about targets
# and other properties processed in source and environments trees before
# these tasks can be performed.
#
#----------------------------------------------------------------------------
# Installation of optional read-only architecture independent data files.
# E.g. Examples, data libraries, documentation.
# Done before toolchain generation because it may affect what we have to do
# there!
#
include(Geant4InstallData)

#------------------------------------------------------------------------------
# Generate any Use/Config files here once everything else has been processed
# e.g. "UseGeant4.cmake", "Geant4Config.cmake" - library dependencies etc
#
include(Geant4ToolchainBackwardCompatibility)
include(Geant4ConfigureConfigScript)

#----------------------------------------------------------------------------
# Create the Geant4Config files and supporting modules
#
include(Geant4BuildProjectConfig)

#------------------------------------------------------------------------------
# Add the examples and tests subdirectories
# This is done after the Geant4Config.cmake file has been generated
if(GEANT4_ENABLE_TESTING)
  include(Geant4CTest)
  add_subdirectory(tests)
  if(EXISTS ${CMAKE_SOURCE_DIR}/benchmarks)
    add_subdirectory(benchmarks)
  endif()
endif()

if(GEANT4_BUILD_EXAMPLES)
  set(Geant4_DIR ${CMAKE_BINARY_DIR} CACHE PATH "Current build directory")
  add_subdirectory(examples)
endif()

# - Install example code to datarootdir
install(DIRECTORY examples
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Geant4-${Geant4_VERSION}
  COMPONENT Examples
  PATTERN "CVS" EXCLUDE
  PATTERN ".svn" EXCLUDE
  )


#-----------------------------------------------------------------------------
# Provide packaging with CPack.
include(Geant4CPackBase)

#----------------------------------------------------------------------------
# Final output - show what's been enabled so that user knows what's
# happening - also useful for later problem solving!
#
GEANT4_PRINT_ENABLED_FEATURES()

