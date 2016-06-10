#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4geometry
# Package: Geant4.src.G4analysis
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 72312 2013-07-15 16:01:02Z gcosmo $
#
#------------------------------------------------------------------------------

add_subdirectory(csv)
add_subdirectory(hntools)
add_subdirectory(management)
add_subdirectory(root)
add_subdirectory(xml)
add_subdirectory(g4tools)

if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4analysis
        COMPONENTS
            csv/sources.cmake
            hntools/sources.cmake
            management/sources.cmake
            root/sources.cmake
            xml/sources.cmake
)
endif()
