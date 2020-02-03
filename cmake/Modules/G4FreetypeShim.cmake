# - G4FreetypeShim
#
# Geant4's Geant4Config.cmake file aims to support CMake 3.8 and newer
# The Freetype dependency is located through CMake's builtin FindFreetype
# module and linked through the Freetype::Freetype imported target.
# This target is however only available from CMake 3.10, so recreate
# Freetype::Freetype target if it does not exist
if(Freetype_FOUND)
  if(NOT TARGET Freetype::Freetype)
    add_library(Freetype::Freetype UNKNOWN IMPORTED)
    set_target_properties(Freetype::Freetype PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${FREETYPE_INCLUDE_DIRS}")

    if(FREETYPE_LIBRARY_RELEASE)
      set_property(TARGET Freetype::Freetype APPEND PROPERTY
        IMPORTED_CONFIGURATIONS RELEASE)
      set_target_properties(Freetype::Freetype PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
        IMPORTED_LOCATION_RELEASE "${FREETYPE_LIBRARY_RELEASE}")
    endif()

    if(FREETYPE_LIBRARY_DEBUG)
      set_property(TARGET Freetype::Freetype APPEND PROPERTY
        IMPORTED_CONFIGURATIONS DEBUG)
      set_target_properties(Freetype::Freetype PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
        IMPORTED_LOCATION_DEBUG "${FREETYPE_LIBRARY_DEBUG}")
    endif()

    if(NOT FREETYPE_LIBRARY_RELEASE AND NOT FREETYPE_LIBRARY_DEBUG)
      set_target_properties(Freetype::Freetype PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${FREETYPE_LIBRARY}")
    endif()
  endif()
endif()

