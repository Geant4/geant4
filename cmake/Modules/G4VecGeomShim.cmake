# Create imported target to help use VecGeom (a mess because the VecGeomConfig is...)
if(VecGeom_FOUND)
  if(NOT TARGET VecGeom::vecgeom)
    add_library(VecGeom::vecgeom UNKNOWN IMPORTED)
    # VecGeom's config file is awful, so the following is neccessary...
    foreach(__vglib ${VECGEOM_LIBRARIES})
      if(__vglib MATCHES ".*libvecgeom\.(a|so|dylib|lib|dll)$")
        set(VECGEOM_LIBRARY "${__vglib}")
      endif()
    endforeach()

    string(REGEX REPLACE "^\-D|;-D" ";" VECGEOM_COMPILE_DEFINITIONS "${VECGEOM_DEFINITIONS}")
    set_target_properties(VecGeom::vecgeom PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "${VECGEOM_COMPILE_DEFINITIONS}"
      INTERFACE_INCLUDE_DIRECTORIES "${VECGEOM_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${VECGEOM_LIBRARIES}"
      IMPORTED_LOCATION "${VECGEOM_LIBRARY}")
  endif()
endif()

