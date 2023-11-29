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
include(G4OptionalComponents)

# - Provide interface to control use of UI/Vis components
#   Written in a separate module from other optional components because
#   there are many complex options to handle.
include(G4InterfaceOptions)

# - Installation of optional read-only architecture independent data files.
include(G4InstallData)

# - Include testing up front so both source and environments can use it if required
include(G4CTest)

#-----------------------------------------------------------------------
# Add the source and environments subdirectories
# source       : Process all the Geant4 core targets
add_subdirectory(source)

if(GEANT4_USE_PYTHON)
  message(WARNING "Geant4Py is no longer distributed with the Geant4 toolkit."
  " It can be downloaded from:\n  https://github.com/koichi-murakami/g4python \n"
  "and requests for support on installation and use should be directed there.\n")
endif()

#-----------------------------------------------------------------------
# - Perform all post build tasks
# At the CMake level, this simply means that we must know about targets
# and other properties processed in source and environments trees before
# these tasks can be performed.
#
# - Generate any Use/Config/Support files here once everything else has
# - Geant4Make.gmk
include(G4ConfigureGNUMakeHelpers)

# - Pkg-Config and geant4-config
include(G4ConfigurePkgConfigHelpers)

# - Geant4Config.cmake et al
include(G4ConfigureCMakeHelpers)

#-----------------------------------------------------------------------
# - Testing configuration.
# Done here, as projects under 'tests' require Geant4Config.
if(GEANT4_ENABLE_TESTING)
  if(EXISTS ${PROJECT_SOURCE_DIR}/tests)
    add_subdirectory(tests)
  endif()
  if(EXISTS ${PROJECT_SOURCE_DIR}/benchmarks)
    add_subdirectory(benchmarks)
  endif()
  if(EXISTS ${PROJECT_SOURCE_DIR}/verification)
    add_subdirectory(verification)
  endif()
endif()

#-----------------------------------------------------------------------
# - Examples build/install
# ON by default for end users. Developers can switch this OFF if they
# need to save time/space
option(GEANT4_INSTALL_EXAMPLES "Install code and documentation for Geant4 examples" ON)
mark_as_advanced(GEANT4_INSTALL_EXAMPLES)

if(GEANT4_INSTALL_EXAMPLES)
  install(DIRECTORY examples
    DESTINATION ${CMAKE_INSTALL_DATADIR}
    COMPONENT Examples
    PATTERN "CVS" EXCLUDE
    PATTERN ".svn" EXCLUDE
    )
endif()

#-----------------------------------------------------------------------
# - CPack-aging
include(G4CPack)

#-----------------------------------------------------------------------
# Final output - show what's been enabled so that user knows what's
# happening - also useful for later problem solving!
#
geant4_print_enabled_features()

