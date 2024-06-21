#[=======================================================================[.rst:
FindPythia6
---------

Find the Pythia6 event generator library.

Imported Targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` targets:

``Pythia6::Pythia6``
  The Pythia6 ``pythia6`` library, if found.

Result Variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project:

``Pythia6_FOUND``
  true if the Pythia6 headers and libraries were found.

Hints
^^^^^

A user may set ``Pythia6_ROOT`` to a Pythia6 installation root to tell this
module where to look.

#]=======================================================================]

# WE DO NOT ALLOW USE OF THIS MODULE OUTSIDE GEANT4 AND EXAMPLES
if(NOT PROJECT_NAME MATCHES "Geant4|HepMCEx0[1-2]|decayer6")
  message(FATAL_ERROR "This FindPythia6.cmake module is only supported for use in Geant4's HepMCEx0{1,2} "
    "and decayer6 extended examples. No support or extension of this module for use outside of these examples "
    "will be provided.")
endif()

find_library(Pythia6_LIBRARY
  NAMES Pythia6 pythia6 pythia6-$ENV{PYTHIA6_VERSION}
  HINTS $ENV{PYTHIA6} $ENV{PYTHIA6}/lib $ENV{PYTHIA6}/lib64 ${Pythia6_ROOT}/lib ${Pythia6_ROOT}/lib64)

# handle the QUIETLY and REQUIRED arguments and set Pythia6_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Pythia6 DEFAULT_MSG Pythia6_LIBRARY)

mark_as_advanced(Pythia6_LIBRARY)

if(Pythia6_FOUND)
  set(Pythia6_LIBRARIES ${Pythia6_LIBRARY})
  if(NOT TARGET Pythia6::Pythia6)
    add_library(Pythia6::Pythia6 UNKNOWN IMPORTED)
    set_target_properties(Pythia6::Pythia6 PROPERTIES
      IMPORTED_LOCATION ${Pythia6_LIBRARY})
  endif()
endif()
