#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_proc_ci
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_chips.G4hadronic_proc_ci
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 18:58:29 bmorgan Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_proc_ci
    HEADERS
        G4QAtomicElectronScattering.hh
        G4QCaptureAtRest.hh
        G4QCoherentChargeExchange.hh
        G4QDiffraction.hh
        G4QDiscProcessMixer.hh
        G4QDiscreteProcessVector.hh
        G4QElastic.hh
        G4QInelastic.hh
        G4QIonIonElastic.hh
        G4QLowEnergy.hh
        G4QPDGToG4Particle.hh
        G4QSynchRad.hh
    SOURCES
        G4QAtomicElectronScattering.cc
        G4QCaptureAtRest.cc
        G4QCoherentChargeExchange.cc
        G4QDiffraction.cc
        G4QDiscProcessMixer.cc
        G4QElastic.cc
        G4QInelastic.cc
        G4QIonIonElastic.cc
        G4QLowEnergy.cc
        G4QPDGToG4Particle.cc
        G4QSynchRad.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4hadronic_body_ci
        G4hadronic_crosec_ci
        G4hadronic_fragm_ci
        G4ions
        G4leptons
        G4magneticfield
        G4materials
        G4mesons
        G4navigation
        G4partman
        G4procman
        G4shortlived
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

