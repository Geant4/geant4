#[=======================================================================[.rst:
FindPythia8
---------

Find the Pythia8 event generator headers and library.

Imported Targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` targets:

``Pythia8::Pythia8``
  The Pythia8 ``pythia8`` library, if found.

Result Variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project:

``Pythia8_FOUND``
  true if the Pythia8 headers and libraries were found.

Hints
^^^^^

A user may set ``Pythia8_ROOT`` to a Pythia8 installation root to tell this
module where to look.

#]=======================================================================]

# WE ONLY ALLOW USE OF THIS MODULE IN GEANT4 OR py8decayer EXAMPLE
if(NOT PROJECT_NAME MATCHES "Geant4|py8decayer")
  message(FATAL_ERROR "This FindPythia8.cmake module is only supported for use in Geant4's py8decayer "
    "extended example. No support or extension of this module for use outside of this example "
    "will be provided.")
endif()

# Look for the header file
find_path(Pythia8_INCLUDE_DIR 
  NAMES Pythia8/Pythia.h
  HINTS ${Pythia8_ROOT}/include)

if(NOT Pythia8_LIBRARY)
  find_library(Pythia8_LIBRARY
    NAMES pythia8 Pythia8 
    HINTS ${Pythia8_ROOT}/lib ${Pythia8_ROOT}/lib64)
endif()

# Determine the version
if(Pythia8_INCLUDE_DIR AND EXISTS "${Pythia8_INCLUDE_DIR}/Pythia8/Pythia.h")
  file(STRINGS "${Pythia8_INCLUDE_DIR}/Pythia8/Pythia.h" PYTHIA8_H REGEX "^#define PYTHIA_VERSION.*$")
  if(PYTHIA8_H MATCHES "PYTHIA_VERSION ([0-9]\\.[0-9]+)$")
    set(Pythia8_VERSION "${CMAKE_MATCH_1}")
  else()
    set(Pythia8_VERSION "")
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set Pythia8_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Pythia8 
  REQUIRED_VARS Pythia8_INCLUDE_DIR Pythia8_LIBRARY
  VERSION_VAR Pythia8_VERSION)

mark_as_advanced(Pythia8_INCLUDE_DIR Pythia8_LIBRARY)

# Create imported target
if(Pythia8_FOUND)
  if(NOT TARGET Pythia8::Pythia8)
    add_library(Pythia8::Pythia8 UNKNOWN IMPORTED)
    set_target_properties(Pythia8::Pythia8 PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${Pythia8_INCLUDE_DIR}"
      IMPORTED_LOCATION ${Pythia8_LIBRARY})
  endif()
endif()
