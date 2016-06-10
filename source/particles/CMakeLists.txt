#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4particles
# Package: Geant4.src.G4particles
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66892 2013-01-17 10:57:59Z gunter $
#
#------------------------------------------------------------------------------

# Add allocation export symbol for the particles category
add_definitions(-DG4PARTICLES_ALLOC_EXPORT)


add_subdirectory(adjoint)
add_subdirectory(bosons)
add_subdirectory(hadrons)
add_subdirectory(leptons)
add_subdirectory(management)
add_subdirectory(shortlived)
add_subdirectory(utils)

if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4particles
        COMPONENTS
            adjoint/sources.cmake
            bosons/sources.cmake
            hadrons/barions/sources.cmake
            hadrons/ions/sources.cmake
            hadrons/mesons/sources.cmake
            leptons/sources.cmake
            management/sources.cmake
            shortlived/sources.cmake
            utils/sources.cmake
)
endif()

