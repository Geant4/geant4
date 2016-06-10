#------------------------------------------------------------------------------
# sources.cmake
# Module : G4shortlived
# Package: Geant4.src.G4particles.G4shortlived
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 66892 2013-01-17 10:57:59Z gunter $
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
GEANT4_DEFINE_MODULE(NAME G4shortlived
    HEADERS
        G4DiQuarks.hh
        G4ExcitedBaryonConstructor.hh
        G4ExcitedBaryons.hh
        G4ExcitedDeltaConstructor.hh
        G4ExcitedLambdaConstructor.hh
        G4ExcitedMesonConstructor.hh
        G4ExcitedMesons.hh
        G4ExcitedNucleonConstructor.hh
        G4ExcitedSigmaConstructor.hh
        G4ExcitedXiConstructor.hh
        G4Gluons.hh
        G4Quarks.hh
        G4ShortLivedConstructor.hh
        G4VShortLivedParticle.hh
    SOURCES
        G4DiQuarks.cc
        G4ExcitedBaryonConstructor.cc
        G4ExcitedBaryons.cc
        G4ExcitedDeltaConstructor.cc
        G4ExcitedLambdaConstructor.cc
        G4ExcitedMesonConstructor.cc
        G4ExcitedMesons.cc
        G4ExcitedNucleonConstructor.cc
        G4ExcitedSigmaConstructor.cc
        G4ExcitedXiConstructor.cc
        G4Gluons.cc
        G4Quarks.cc
        G4ShortLivedConstructor.cc
        G4VShortLivedParticle.cc
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

