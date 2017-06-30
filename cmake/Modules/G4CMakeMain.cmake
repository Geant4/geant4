#.rst:
# G4CMakeMain
# -----------
#
# Main script for Geant4 CMake build scripting.
#
# Stored separately from top level script purely for clarity in tagging.
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

# Include once and fail if already included
if(NOT __G4CMAKEMAIN_INCLUDED)
  set(__G4CMAKEMAIN_INCLUDED TRUE)
else()
  message(FATAL_ERROR "G4CMakeMain can only be included once")
endif()

# Shouldn't need to set CMAKE_MODULE_PATH here as done at top level.
# Review whether this is a more appropriate location.

#-----------------------------------------------------------------------
# - Add functionality provided by standard and custom modules
# See the documentation in each of these modules for further details.

# - Core CMake settings and options
include(G4CMakeSettings)

# - Core Compile/Link settings and options
include(G4BuildSettings)

# - Core API
include(G4DeveloperAPI)

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
# - Geant4Make.gmk
include(G4ConfigureGNUMakeHelpers)

# - Pkg-Config and geant4-config
include(G4ConfigurePkgConfigHelpers)

# - Geant4Config.cmake et al
include(G4ConfigureCMakeHelpers)

#-----------------------------------------------------------------------
# - Testing configuration.
# Done here, as projects under 'tests' require Geant4Config.
include(Geant4CTest)
if(GEANT4_ENABLE_TESTING)
  add_subdirectory(tests)
  if(EXISTS ${CMAKE_SOURCE_DIR}/benchmarks)
    add_subdirectory(benchmarks)
  endif()
  if(EXISTS ${CMAKE_SOURCE_DIR}/verification)
    add_subdirectory(verification)
  endif()
endif()

#-----------------------------------------------------------------------
# - Examples build/install
# NB: Build of examples is a *testing* proceedure. It is *not* intended
# that examples be built and installed as part of a full Geant4 install.
option(GEANT4_BUILD_EXAMPLES "Build all the examples of the project" OFF)
GEANT4_ADD_FEATURE(GEANT4_BUILD_EXAMPLES "Build all the examples of the project")
mark_as_advanced(GEANT4_BUILD_EXAMPLES)

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

