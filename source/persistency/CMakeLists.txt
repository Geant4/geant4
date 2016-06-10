#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4persistency
# Package: Geant4.src.G4persistency
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66892 2013-01-17 10:57:59Z gunter $
#
#------------------------------------------------------------------------------

add_subdirectory(ascii)

# GDML only if required
if(GEANT4_USE_GDML)
    add_subdirectory(gdml)
endif()

add_subdirectory(mctruth)


if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    if(GEANT4_USE_GDML)
        GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4persistency
            COMPONENTS
                ascii/sources.cmake
                gdml/sources.cmake
                mctruth/sources.cmake
        )
    else()
        GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4persistency
            COMPONENTS
                ascii/sources.cmake
                mctruth/sources.cmake
        )
    endif()
endif()

