#-----------------------------------------------------------------------
# - Top Level CMakeLists.txt for Geant4 Build
#
# 21st September 2010 Ben Morgan
#
# $Id: CMakeLists.txt 106587 2017-10-13 13:07:50Z gcosmo $
#

#-----------------------------------------------------------------------
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

#-----------------------------------------------------------------------
# - Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

# - Any policy requirements should go here

set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
   ${CMAKE_SOURCE_DIR}/cmake/Modules/Geant4MakeRules_cxx.cmake)

#-----------------------------------------------------------------------
# - Project definition and basic configuration
#
project(Geant4)

# - Versioning. We do this here for now
set(${PROJECT_NAME}_VERSION "10.3.3")
set(${PROJECT_NAME}_VERSION_MAJOR "10")
set(${PROJECT_NAME}_VERSION_MINOR "3")
set(${PROJECT_NAME}_VERSION_PATCH "3")


# - Prepend our own CMake Modules to the search path
# NB: if our custom modules include others that we don't supply, those in
# the base path will be used, so watch for incompatibilities!!
#
set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake/Modules
    ${CMAKE_MODULE_PATH})

#-----------------------------------------------------------------------
# - Add functionality provided by standard and custom modules
# See the documentation in each of these modules for further details.
#
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

#-----------------------------------------------------------------------
# Add the source and environments subdirectories
# source       : Process all the Geant4 core targets
# environments : Process optional wrappings of Geant4 (NOTYETIMPLEMENTED)
add_subdirectory(source)
#add_subdirectory(environments)

#-----------------------------------------------------------------------
# - Perform all post build tasks
# At the CMake level, this simply means that we must know about targets
# and other properties processed in source and environments trees before
# these tasks can be performed.
#
# - Installation of optional read-only architecture independent data files.
# E.g. Examples, data libraries, documentation.
# Done before toolchain generation because it may affect what we have to do
# there!
#
include(Geant4InstallData)

# - Generate any Use/Config/Support files here once everything else has
# been processed e.g. "UseGeant4.cmake", "Geant4Config.cmake", library
# dependencies etc.
# - Geant4Make
include(Geant4ToolchainBackwardCompatibility)

# - 'geant4-config'
include(Geant4ConfigureConfigScript)

# - Geant4Config.cmake
include(Geant4BuildProjectConfig)

#-----------------------------------------------------------------------
# - Testing configuration.
# Done here, as projects under 'tests' require Geant4Config.
if(GEANT4_ENABLE_TESTING)
  include(Geant4CTest)
  add_subdirectory(tests)
  if(EXISTS ${CMAKE_SOURCE_DIR}/benchmarks)
    add_subdirectory(benchmarks)
  endif()
endif()

#-----------------------------------------------------------------------
# - Examples build/install
# NB: Build of examples is a *testing* proceedure. It is *not* intended
# that examples be built and installed as part of a full Geant4 install.
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

#-----------------------------------------------------------------------
# - CPack-aging
include(Geant4CPackBase)

#-----------------------------------------------------------------------
# Final output - show what's been enabled so that user knows what's
# happening - also useful for later problem solving!
#
GEANT4_PRINT_ENABLED_FEATURES()

