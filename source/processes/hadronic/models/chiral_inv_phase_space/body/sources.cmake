#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_body_ci
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_chips.G4hadronic_body_ci
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 18:57:25 bmorgan Exp $
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

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_body_ci
    HEADERS
        G4Chips.hh
        G4QBesIKJY.hh
        G4QCHIPSWorld.hh
        G4QCandidate.hh
        G4QCandidateVector.hh
        G4QChipolino.hh
        G4QContent.hh
        G4QDecayChan.hh
        G4QDecayChanVector.hh
        G4QEnvironment.hh
        G4QException.hh
        G4QHadron.hh
        G4QHadronVector.hh
        G4QInteraction.hh
        G4QInteractionVector.hh
        G4QIsotope.hh
        G4QNucleus.hh
        G4QPDGCode.hh
        G4QPDGCodeVector.hh
        G4QParentCluster.hh
        G4QParentClusterVector.hh
        G4QParticle.hh
        G4QParticleVector.hh
        G4QParton.hh
        G4QPartonPair.hh
        G4QPartonPairVector.hh
        G4QPartonVector.hh
        G4QProbability.hh
        G4QString.hh
        G4QStringVector.hh
        G4Quasmon.hh
        G4QuasmonVector.hh
    SOURCES
        G4QBesIKJY.cc
        G4QCHIPSWorld.cc
        G4QCandidate.cc
        G4QChipolino.cc
        G4QContent.cc
        G4QDecayChan.cc
        G4QEnvironment.cc
        G4QException.cc
        G4QHadron.cc
        G4QInteraction.cc
        G4QIsotope.cc
        G4QNucleus.cc
        G4QPDGCode.cc
        G4QParentCluster.cc
        G4QParticle.cc
        G4QParton.cc
        G4QPartonPair.cc
        G4QProbability.cc
        G4QString.cc
        G4Quasmon.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4globman
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

