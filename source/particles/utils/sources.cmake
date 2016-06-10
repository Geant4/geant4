#------------------------------------------------------------------------------
# sources.cmake
# Module : G4partutils
# Package: Geant4.src.G4particles.G4partutils
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
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4partutils
    HEADERS
        G4HtmlPPReporter.hh
        G4IsotopeMagneticMomentTable.hh
        G4SimplePPReporter.hh
        G4TextPPReporter.hh
        G4TextPPRetriever.hh
        G4VParticlePropertyReporter.hh
        G4VParticlePropertyRetriever.hh
    SOURCES
        G4HtmlPPReporter.cc
        G4IsotopeMagneticMomentTable.cc
        G4SimplePPReporter.cc
        G4TextPPReporter.cc
        G4TextPPRetriever.cc
        G4VParticlePropertyReporter.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4intercoms
        G4materials
        G4partman
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4intercoms
        G4materials
    LINK_LIBRARIES
)

# List any source specific properties here

