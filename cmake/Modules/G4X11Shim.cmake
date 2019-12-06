# - G4X11Shim
#
# Geant4's Geant4Config.cmake file aims to support CMake 3.8 and newer
# The X11 dependency is located through CMake's builtin FindX11
# module and linked through the X11:: imported targets.
# These targets are however only available from CMake 3.14, so recreate
# those we need if they do not exist.

if(X11_FOUND)
  if (NOT TARGET X11::X11)
      add_library(X11::X11 UNKNOWN IMPORTED)
      set_target_properties(X11::X11 PROPERTIES
        IMPORTED_LOCATION "${X11_X11_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${X11_X11_INCLUDE_PATH}")
  endif ()

  if (X11_ICE_FOUND AND NOT TARGET X11::ICE)
    add_library(X11::ICE UNKNOWN IMPORTED)
    set_target_properties(X11::ICE PROPERTIES
      IMPORTED_LOCATION "${X11_ICE_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${X11_ICE_INCLUDE_PATH}")
  endif ()

  if (X11_SM_FOUND AND NOT TARGET X11::SM)
    add_library(X11::SM UNKNOWN IMPORTED)
    set_target_properties(X11::SM PROPERTIES
      IMPORTED_LOCATION "${X11_SM_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${X11_SM_INCLUDE_PATH}")
  endif ()

  if (X11_Xext_FOUND AND NOT TARGET X11::Xext)
    # CMake < 3.14 won't search for Xext headers
    if(NOT X11_Xext_INCLUDE_PATH)
      find_path(X11_Xext_INCLUDE_PATH X11/extensions/Xext.h ${X11_INC_SEARCH_PATH})
      mark_as_advanced(X11_Xext_INCLUDE_PATH)
    endif()
    add_library(X11::Xext UNKNOWN IMPORTED)
    set_target_properties(X11::Xext PROPERTIES
      IMPORTED_LOCATION "${X11_Xext_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${X11_Xext_INCLUDE_PATH}"
      INTERFACE_LINK_LIBRARIES "X11::X11")
  endif ()

  if (X11_Xmu_FOUND AND NOT TARGET X11::Xmu)
    add_library(X11::Xmu UNKNOWN IMPORTED)
    set_target_properties(X11::Xmu PROPERTIES
      IMPORTED_LOCATION "${X11_Xmu_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${X11_Xmu_INCLUDE_PATH}"
      INTERFACE_LINK_LIBRARIES "X11::Xt;X11::Xext;X11::X11")
  endif ()

  if (X11_Xpm_FOUND AND NOT TARGET X11::Xpm)
    add_library(X11::Xpm UNKNOWN IMPORTED)
    set_target_properties(X11::Xpm PROPERTIES
      IMPORTED_LOCATION "${X11_Xpm_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${X11_Xpm_INCLUDE_PATH}"
      INTERFACE_LINK_LIBRARIES "X11::X11")
  endif ()

  if (X11_Xt_FOUND AND NOT TARGET X11::Xt)
    add_library(X11::Xt UNKNOWN IMPORTED)
    set_target_properties(X11::Xt PROPERTIES
      IMPORTED_LOCATION "${X11_Xt_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${X11_Xt_INCLUDE_PATH}"
      INTERFACE_LINK_LIBRARIES "X11::ICE;X11::SM;X11::X11")
  endif ()
endif()
