# - G4MotifShim
#
# Geant4's Geant4Config.cmake file aims to support CMake 3.8 and newer
# The Motif dependency is located through CMake's builtin FindMotif
# This module does not support imported targets as of CMake 3.16, so
# provide a target Motif::Xm
if(Motif_FOUND)
  if(NOT TARGET Motif::Xm)
    add_library(Motif::Xm UNKNOWN IMPORTED)
    set_target_properties(Motif::Xm PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${MOTIF_INCLUDE_DIR}"
      IMPORTED_LOCATION "${MOTIF_LIBRARIES}"
    )
  endif()
endif()