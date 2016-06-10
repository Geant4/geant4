#------------------------------------------------------------------------------
# sources.cmake
# Module : G4bosons
# Package: Geant4.src.G4particles.G4bosons
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 75122 2013-10-28 09:51:40Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4bosons
    HEADERS
        G4BosonConstructor.hh
        G4ChargedGeantino.hh
        G4Gamma.hh
        G4Geantino.hh
        G4OpticalPhoton.hh
	G4PhononLong.hh
	G4PhononTransFast.hh
	G4PhononTransSlow.hh
        G4UnknownParticle.hh
    SOURCES
        G4BosonConstructor.cc
        G4ChargedGeantino.cc
        G4Gamma.cc
        G4Geantino.cc
        G4OpticalPhoton.cc
	G4PhononLong.cc
	G4PhononTransFast.cc
	G4PhononTransSlow.cc
        G4UnknownParticle.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4materials
        G4partman
    GLOBAL_DEPENDENCIES
        G4global
        G4materials
    LINK_LIBRARIES
)

# List any source specific properties here

