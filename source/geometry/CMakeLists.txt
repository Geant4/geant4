#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4geometry
# Package: Geant4.src.G4geometry
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 75578 2013-11-04 12:03:33Z gcosmo $
#
#------------------------------------------------------------------------------

# Add allocation export symbol for the geometry category
add_definitions(-DG4GEOM_ALLOC_EXPORT)

add_subdirectory(biasing)
add_subdirectory(divisions)
add_subdirectory(magneticfield)
add_subdirectory(management)
add_subdirectory(navigation)
add_subdirectory(solids)
add_subdirectory(volumes)

if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4geometry
        COMPONENTS
            biasing/sources.cmake
            divisions/sources.cmake
            magneticfield/sources.cmake
            management/sources.cmake
            navigation/sources.cmake
            solids/Boolean/sources.cmake
            solids/CSG/sources.cmake
            solids/specific/sources.cmake
            solids/usolids/sources.cmake
            volumes/sources.cmake
)
endif()

