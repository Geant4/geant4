#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4globman
# Package: Geant4.src.G4global.G4globman
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
#
#------------------------------------------------------------------------------

execute_process(
    COMMAND git rev-parse HEAD
    OUTPUT_VARIABLE GIT_SHA1
)

string(REGEX REPLACE "\n$" "" GIT_SHA1 "${GIT_SHA1}")

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4Version.hh
    ${CMAKE_CURRENT_SOURCE_DIR}/include/G4Version.hh
)

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

