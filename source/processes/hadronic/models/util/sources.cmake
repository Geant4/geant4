#------------------------------------------------------------------------------
# sources.cmake
# Module : G4had_mod_util
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_mod_util
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 104985 2017-07-03 15:14:26Z gcosmo $
#
# 20110727  M. Kelsey -- Add G4DecayKineticTracks
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4had_mod_util
    HEADERS
        G4Clebsch.hh
	G4DecayKineticTracks.hh
        G4DecayStrongResonances.hh
        G4ExcitedString.hh
        G4ExcitedStringVector.hh
        G4Fancy3DNucleus.hh
        G4Fancy3DNucleusHelper.hh
        G4FermiMomentum.hh
        G4Fragment.hh
        G4FragmentVector.hh
        G4GeneralPhaseSpaceDecay.hh
	G4HadDecayGenerator.hh
	G4HadPhaseSpaceGenbod.hh
	G4HadPhaseSpaceKopylov.hh
	G4HadPhaseSpaceNBodyAsai.hh
        G4KineticTrack.hh
        G4KineticTrackVector.hh
        G4LegendrePolynomial.hh
        G4NuclearFermiDensity.hh
        G4NuclearPolarizationStore.hh
        G4NuclearShellModelDensity.hh
        G4Nucleon.hh
        G4Parton.hh
        G4PartonVector.hh
        G4PolynomialPDF.hh
        G4SampleResonance.hh
	G4VHadDecayAlgorithm.hh
	G4VHadPhaseSpaceAlgorithm.hh
        G4WilsonRadius.hh
    SOURCES
        G4Clebsch.cc
	G4DecayKineticTracks.cc
        G4DecayStrongResonances.cc
        G4ExcitedString.cc
        G4Fancy3DNucleus.cc
        G4FermiMomentum.cc
        G4Fragment.cc
        G4GeneralPhaseSpaceDecay.cc
	G4HadDecayGenerator.cc
	G4HadPhaseSpaceGenbod.cc
	G4HadPhaseSpaceKopylov.cc
	G4HadPhaseSpaceNBodyAsai.cc
        G4KineticTrack.cc
        G4KineticTrackVector.cc
        G4LegendrePolynomial.cc
        G4NuclearFermiDensity.cc
        G4NuclearPolarizationStore.cc
        G4NuclearShellModelDensity.cc
        G4Nucleon.cc
        G4Parton.cc
        G4PolynomialPDF.cc
        G4SampleResonance.cc
	G4VHadDecayAlgorithm.cc
	G4VHadPhaseSpaceAlgorithm.cc
        G4WilsonRadius.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4hadronic_util
        G4hepnumerics
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4navigation
        G4partman
        G4procman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

