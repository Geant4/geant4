#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_crosec_ci
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_chips.G4hadronic_crosec_ci
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 18:57:41 bmorgan Exp $
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

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_crosec_ci
    HEADERS
        G4QANuANuNuclearCrossSection.hh
        G4QANuENuclearCrossSection.hh
        G4QANuMuNuclearCrossSection.hh
        G4QAntiBaryonElasticCrossSection.hh
        G4QAntiBaryonNuclearCrossSection.hh
        G4QAntiBaryonPlusNuclearCrossSection.hh
        G4QDiffractionRatio.hh
        G4QElectronNuclearCrossSection.hh
        G4QHyperonElasticCrossSection.hh
        G4QHyperonNuclearCrossSection.hh
        G4QHyperonPlusElasticCrossSection.hh
        G4QHyperonPlusNuclearCrossSection.hh
        G4QIonIonCrossSection.hh
        G4QKaonMinusElasticCrossSection.hh
        G4QKaonMinusNuclearCrossSection.hh
        G4QKaonPlusElasticCrossSection.hh
        G4QKaonPlusNuclearCrossSection.hh
        G4QKaonZeroNuclearCrossSection.hh
        G4QMuonNuclearCrossSection.hh
        G4QNeutronCaptureRatio.hh
        G4QNeutronElasticCrossSection.hh
        G4QNeutronNuclearCrossSection.hh
        G4QNuENuclearCrossSection.hh
        G4QNuMuNuclearCrossSection.hh
        G4QNuNuNuclearCrossSection.hh
        G4QPhotonNuclearCrossSection.hh
        G4QPionMinusElasticCrossSection.hh
        G4QPionMinusNuclearCrossSection.hh
        G4QPionPlusElasticCrossSection.hh
        G4QPionPlusNuclearCrossSection.hh
        G4QProtonElasticCrossSection.hh
        G4QProtonNuclearCrossSection.hh
        G4QTauNuclearCrossSection.hh
        G4QuasiFreeRatios.hh
        G4VQCrossSection.hh
    SOURCES
        G4QANuANuNuclearCrossSection.cc
        G4QANuENuclearCrossSection.cc
        G4QANuMuNuclearCrossSection.cc
        G4QAntiBaryonElasticCrossSection.cc
        G4QAntiBaryonNuclearCrossSection.cc
        G4QAntiBaryonPlusNuclearCrossSection.cc
        G4QDiffractionRatio.cc
        G4QElectronNuclearCrossSection.cc
        G4QHyperonElasticCrossSection.cc
        G4QHyperonNuclearCrossSection.cc
        G4QHyperonPlusElasticCrossSection.cc
        G4QHyperonPlusNuclearCrossSection.cc
        G4QIonIonCrossSection.cc
        G4QKaonMinusElasticCrossSection.cc
        G4QKaonMinusNuclearCrossSection.cc
        G4QKaonPlusElasticCrossSection.cc
        G4QKaonPlusNuclearCrossSection.cc
        G4QKaonZeroNuclearCrossSection.cc
        G4QMuonNuclearCrossSection.cc
        G4QNeutronCaptureRatio.cc
        G4QNeutronElasticCrossSection.cc
        G4QNeutronNuclearCrossSection.cc
        G4QNuENuclearCrossSection.cc
        G4QNuMuNuclearCrossSection.cc
        G4QNuNuNuclearCrossSection.cc
        G4QPhotonNuclearCrossSection.cc
        G4QPionMinusElasticCrossSection.cc
        G4QPionMinusNuclearCrossSection.cc
        G4QPionPlusElasticCrossSection.cc
        G4QPionPlusNuclearCrossSection.cc
        G4QProtonElasticCrossSection.cc
        G4QProtonNuclearCrossSection.cc
        G4QTauNuclearCrossSection.cc
        G4QuasiFreeRatios.cc
        G4VQCrossSection.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4globman
        G4hadronic_body_ci
        G4hadronic_crosec_ci
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

