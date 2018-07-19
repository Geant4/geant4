#.rst:
# G4CMakeSettings.cmake
# ---------------------
#
# Set defaults for core CMake behaviour useful or required for building,
# testing and installing Geant4.

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

if(NOT __G4CMAKESETTINGS_INCLUDED)
  set(__G4CMAKESETTINGS_INCLUDED TRUE)
else()
  return()
endif()

#-----------------------------------------------------------------------
#.rst:
# General Configuration Settings
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# The following CMake modules and variables are used to control configuration time
# behaviour such as location of dependencies.
#

#.rst:
# - ``G4CMakeUtilities`
#
#   - Included to provide additional macros and functions
#
include(G4CMakeUtilities)

#
# - ``CMAKE_EXPORT_NO_PACKAGE_REGISTRY`` : ON
# - ``CMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY`` : ON
# - ``CMAKE_FIND_PACKAGE_NO_SYSTEM_PACKAGE_REGISTRY`` : ON
#
#   - These variables are set to ensure that CMake never creates or searches
#     for config files in any package registry. This prevents the
#     :cmake:command:`find_package <cmake:command:find_package>`
#     command from locating potentially spurious config files.
#
set(CMAKE_EXPORT_NO_PACKAGE_REGISTRY ON)
set(CMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY ON)
set(CMAKE_FIND_PACKAGE_NO_SYSTEM_PACKAGE_REGISTRY ON)

#-----------------------------------------------------------------------
#.rst:
# General Build Settings
# ^^^^^^^^^^^^^^^^^^^^^^
#
# The following CMake variables and options are configured by default
# when including this module:
#
# - ``CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE`` : ON
#
#   - Force project directories to appear first in any list of include paths.
#     This applies to both full paths and those created by generator expressions.
#
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)

#.rst:
# - ``CMAKE_LINK_DEPENDS_NO_SHARED`` : ON
#
#   - Do not relink a target to any shared library dependencies when
#     only the shared library implementation has changed.
#
set(CMAKE_LINK_DEPENDS_NO_SHARED ON)

#.rst:
# A custom ``validate_sources`` targets is declared to compare files listed
# in `sources.cmake` scripts with on-disk files. As sources.cmake lists
# source files of Geant4 explicitely, we can often a mismatch can occur
# between this list and what's actually on disk.
#
# The ``validate_sources`` target executes a CMake script to
# check for these errors and report mismatches. It fails with FATAL_ERROR
# if any mismatch is found, but will not do so until it has reported
# all errors.
# Configure the script
configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/Templates/geant4_validate_sources.cmake.in
  ${PROJECT_BINARY_DIR}/geant4_validate_sources.cmake
  @ONLY
  )

# Create the target
add_custom_target(validate_sources
  COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BINARY_DIR}/geant4_validate_sources.cmake
  COMMENT "Validating Geant4 Module Source Lists"
  )

#-----------------------------------------------------------------------
#.rst:
# General Installation Settings
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# CMake's builtin `GNUInstallDirs` module is used to set and provide variables
# for the destinations to which for executables, libraries and other files
# should be installed.
#
include(GNUInstallDirs)

# Check for non-relocatable install directories
foreach(dir
    BINDIR
    LIBDIR
    INCLUDEDIR
    DATAROOTDIR
    DATADIR
    MANDIR
    DOCDIR
    )
  if(IS_ABSOLUTE ${CMAKE_INSTALL_${dir}})
    set(CMAKE_INSTALL_IS_NONRELOCATABLE 1)
  endif()
endforeach()

#.rst:
# The following CMake variables are additionally set to control install
# policy and reporting:

#.rst:
# - ``CMAKE_INSTALL_MESSAGE`` : ``LAZY``
#
#   - Only report new or updated files installed by the ``install`` target.
#
set(CMAKE_INSTALL_MESSAGE LAZY)

#.rst:
# An `uninstall` target is provided to assist in removing previously installed
# files. It will not remove any installed directories per standard behaviour.
# Note that if files are removed from the install manifest between invocations
# of the `install` and `uninstall` targets, they will *not* be removed by the
# latter.
#
include(CMakeUninstallTarget)

