#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_deex_fermi_breakup
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_deex.G4hadronic_deex_fermi_breakup
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 18:59:42 bmorgan Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
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
include_directories(${CMAKE_SOURCE_DIR}/source/particles/shortlived/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/cross_sections/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_deex_fermi_breakup
    HEADERS
        G4B9FermiFragment.hh
        G4Be8FermiFragment.hh
        G4FermiBreakUp.hh
        G4FermiConfiguration.hh
        G4FermiConfigurationList.hh
        G4FermiFragmentsPool.hh
        G4FermiPhaseSpaceDecay.hh
        G4He5FermiFragment.hh
        G4Li5FermiFragment.hh
        G4StableFermiFragment.hh
        G4UnstableFermiFragment.hh
        G4VFermiBreakUp.hh
        G4VFermiFragment.hh
    SOURCES
        G4B9FermiFragment.cc
        G4Be8FermiFragment.cc
        G4FermiBreakUp.cc
        G4FermiConfiguration.cc
        G4FermiConfigurationList.cc
        G4FermiFragmentsPool.cc
        G4FermiPhaseSpaceDecay.cc
        G4He5FermiFragment.cc
        G4Li5FermiFragment.cc
        G4StableFermiFragment.cc
        G4UnstableFermiFragment.cc
        G4VFermiBreakUp.cc
        G4VFermiFragment.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4globman
        G4had_mod_util
        G4hadronic_mgt
        G4hadronic_util
        G4hadronic_xsect
        G4hepnumerics
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4procman
        G4shortlived
        G4track
    GLOBAL_DEPENDENCIES
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

