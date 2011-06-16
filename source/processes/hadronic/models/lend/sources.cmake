#------------------------------------------------------------------------------
# sources.cmake
# Module : G4had_lend
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_lend
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-11-23 13:58:03 gcosmo Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${EXPAT_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/cross_sections/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/photon_evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/low_energy/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/processes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)


#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4had_lend
    HEADERS
        G4GIDI.hh
        G4GIDI_Misc.hh
        G4GIDI_map.hh
        G4GIDI_mass.hh
        G4GIDI_target.hh
        G4LENDCapture.hh
        G4LENDCaptureCrossSection.hh
        G4LENDCrossSection.hh
        G4LENDElastic.hh
        G4LENDElasticCrossSection.hh
        G4LENDFission.hh
        G4LENDFissionCrossSection.hh
        G4LENDHeader.hh
        G4LENDInelastic.hh
        G4LENDInelasticCrossSection.hh
        G4LENDManager.hh
        G4LENDModel.hh
        G4LENDUsedTarget.hh
        gString.h
        statusMessageReporting.h
        tpi_IDs.h
        tpia_defs.h
        tpia_depot.h
        tpia_map.h
        tpia_mass.h
        tpia_misc.h
        tpia_target.h
        xData.h
        xDataExtras.h
    SOURCES
        G4GIDI.cc
        G4GIDI_Misc.cc
        G4GIDI_map.cc
        G4GIDI_mass.cc
        G4GIDI_target.cc
        G4LENDCapture.cc
        G4LENDCaptureCrossSection.cc
        G4LENDCrossSection.cc
        G4LENDElastic.cc
        G4LENDElasticCrossSection.cc
        G4LENDFission.cc
        G4LENDFissionCrossSection.cc
        G4LENDInelastic.cc
        G4LENDInelasticCrossSection.cc
        G4LENDManager.cc
        G4LENDModel.cc
        G4LENDUsedTarget.cc
        gString.cc
        statusMessageReporting.cc
        tpi_IDs.cc
        tpia_angular.cc
        tpia_angularEnergy.cc
        tpia_channel.cc
        tpia_decayChannel.cc
        tpia_depot.cc
        tpia_frame.cc
        tpia_kinetics.cc
        tpia_Legendre.cc
        tpia_map.cc
        tpia_mass.cc
        tpia_misc.cc
        tpia_multiplicity.cc
        tpia_particle.cc
        tpia_product.cc
        tpia_samplingMethods.cc
        tpia_target.cc
        tpia_target_heated.cc
        xData.cc
        xDataExtras.cc
        xDataMisc.cc
        xData_1d_x.cc
        xData_2d_xindex_y.cc
        xData_2d_xshared_yhistogram.cc
        xData_2d_xy.cc
        xData_matrix.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4hadronic_LE
        G4hadronic_deex_management
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util
        G4hadronic_mgt
        G4hadronic_proc
        G4hadronic_util
        G4hadronic_xsect
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
        ${EXPAT_LIBRARIES}
)

# List any source specific properties here
#
# Not this again... Need to add *another* compile definition to lend sources
# for DLLs....
# We do this quick and dirty, add the definition directly. I *don't want it 
# polluting compile of all other processes modules
# in global mode
#
set_source_files_properties(
    ${G4had_lend_SOURCES}
    PROPERTIES COMPILE_DEFINITIONS G4PROCESSES_EXPORT
)


