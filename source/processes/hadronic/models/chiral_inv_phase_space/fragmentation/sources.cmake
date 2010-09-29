#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_fragm_ci
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_chips.G4hadronic_fragm_ci
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 18:57:57 bmorgan Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/shortlived/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/chiral_inv_phase_space/body/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/chiral_inv_phase_space/cross_sections/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/chiral_inv_phase_space/fragmentation/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_fragm_ci
    HEADERS
        G4QFragmentation.hh
        G4QIonIonCollision.hh
    SOURCES
        G4QFragmentation.cc
        G4QIonIonCollision.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4globman
        G4hadronic_body_ci
        G4hadronic_crosec_ci
        G4hadronic_fragm_ci
        G4ions
        G4leptons
        G4mesons
        G4partman
        G4shortlived
    GLOBAL_DEPENDENCIES
        G4global
        G4particles
    LINK_LIBRARIES
)

# List any source specific properties here

