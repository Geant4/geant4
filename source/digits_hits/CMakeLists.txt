#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4digits_hits
# Package: Geant4.src.G4digits_hits
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66892 2013-01-17 10:57:59Z gunter $
#
#------------------------------------------------------------------------------

# Add allocation export symbol for the digits_hits category
add_definitions(-DG4DIGI_ALLOC_EXPORT)


add_subdirectory(detector)
add_subdirectory(digits)
add_subdirectory(hits)
add_subdirectory(scorer)
add_subdirectory(utils)

if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4digits_hits
        COMPONENTS
            detector/sources.cmake
            digits/sources.cmake
            hits/sources.cmake
            scorer/sources.cmake
            utils/sources.cmake
)
endif()

