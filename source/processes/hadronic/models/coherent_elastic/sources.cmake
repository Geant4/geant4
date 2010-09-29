#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_coherent_elastic
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_coherent_elastic
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.2 2010-09-29 20:26:23 bmorgan Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/cross_sections/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/chiral_inv_phase_space/body/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/chiral_inv_phase_space/cross_sections/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/neutron_hp/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_coherent_elastic
    HEADERS
        G4CHIPSElastic.hh
        G4CHIPSElasticXS.hh
        G4ChargeExchange.hh
        G4ChargeExchangeProcess.hh
        G4DiffuseElastic.hh
        G4ElasticHadrNucleusHE.hh
        G4HadronElastic.hh
        G4LEnp.hh
        G4LEnpData.hh
        G4LEpp.hh
        G4LEppData.hh
        G4NuclNuclDiffuseElastic.hh
        G4UHadronElasticProcess.hh
        G4VHadronElastic.hh
        G4WHadronElasticProcess.hh
    SOURCES
        G4CHIPSElastic.cc
        G4CHIPSElasticXS.cc
        G4ChargeExchange.cc
        G4ChargeExchangeProcess.cc
        G4DiffuseElastic.cc
        G4ElasticHadrNucleusHE.cc
        G4HadronElastic.cc
        G4LEnp.cc
        G4LEpp.cc
        G4NuclNuclDiffuseElastic.cc
        G4UHadronElasticProcess.cc
        G4VHadronElastic.cc
        G4WHadronElasticProcess.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_neu_hp
        G4hadronic_body_ci
        G4hadronic_crosec_ci
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

