#------------------------------------------------------------------------------
# sources.cmake
# Module : G4had_fission
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_fission
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 08/11/2013
#
# $Id:$
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/particle_hp/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/processes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4had_fission
    HEADERS
        G4FissLib.hh
        G4FissionLibrary.hh
        G4LFission.hh
        G4LLNLFission.hh
        G4fissionEvent.hh
    SOURCES
        G4FissLib.cc
        G4FissionLibrary.cc
        G4LFission.cc
        G4LLNLFission.cc
        G4SmpGEng.cc
        G4SmpIsoDir.cc
        G4SmpNEngCf252.cc
        G4SmpNVel.cc
        G4SmpNuDistDataPu239.cc
        G4SmpNuDistDataPu239_241.cc
        G4SmpNuDistDataPu239_241_MC.cc
        G4SmpNuDistDataU232_234_236_238.cc
        G4SmpNuDistDataU232_234_236_238_MC.cc
        G4SmpNuDistDataU233_235.cc
        G4SmpNuDistDataU233_235_MC.cc
        G4SmpNuDistDataU235.cc
        G4SmpNuDistDataU238.cc
        G4SmpNugDist.cc
        G4SmpPVel.cc
        G4SmpSpNuDistData.cc
        G4SmpSpNubarData.cc
        G4SmpSpNugDistData.cc
        G4SmpTerrell.cc
        G4SmpWatt.cc
        G4fissionEvent.cc
        G4fissionerr.cc
        G4rngc.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_neu_hp
        G4had_par_hp
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
)

# List any source specific properties here

