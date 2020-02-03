# - G4EXPATShim
# 
# Geant4's Geant4Config.cmake file aims to support CMake 3.8 and newer
# The EXPAT dependency is located through CMake's builtin FindEXPAT
# module and linked through the EXPAT::EXPAT imported target.
# This target is however only available from CMake 3.10, so recreate
# EXPAT::EXPAT target if it does not exist
if(EXPAT_FOUND)
  if(NOT TARGET EXPAT::EXPAT)
    add_library(EXPAT::EXPAT UNKNOWN IMPORTED)
    set_target_properties(EXPAT::EXPAT PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      IMPORTED_LOCATION "${EXPAT_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${EXPAT_INCLUDE_DIRS}")
  endif()
endif()