#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4physicslists
# Package: Geant4.src.G4physicslists
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 67969 2013-03-13 09:44:42Z gcosmo $
#
#------------------------------------------------------------------------------

add_subdirectory(builders)
add_subdirectory(lists)
add_subdirectory(constructors)
add_subdirectory(util)

if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4physicslists
        COMPONENTS
            builders/sources.cmake
	    constructors/decay/sources.cmake
	    constructors/electromagnetic/sources.cmake
	    constructors/factory/sources.cmake
	    constructors/gamma_lepto_nuclear/sources.cmake
	    constructors/hadron_elastic/sources.cmake
	    constructors/hadron_inelastic/sources.cmake
	    constructors/ions/sources.cmake
	    constructors/limiters/sources.cmake
	    constructors/stopping/sources.cmake
            lists/sources.cmake
	    util/sources.cmake
)
endif()

