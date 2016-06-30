#------------------------------------------------------------------------------
# sources.cmake
# Module : G4xrays
# Package: Geant4.src.G4processes.G4electromagnetic.G4xrays
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 95093 2016-01-25 08:48:23Z gcosmo $
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/standard/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4xrays
    HEADERS
        G4Cerenkov.hh
        G4ForwardXrayTR.hh
        G4GammaXTRadiator.hh
        G4RegularXTRadiator.hh
        G4Scintillation.hh
        G4ScintillationTrackInformation.hh
        G4StrawTubeXTRadiator.hh
        G4SynchrotronRadiation.hh
        G4SynchrotronRadiationInMat.hh
        G4TransitionRadiation.hh
        G4TransparentRegXTRadiator.hh
        G4VTRModel.hh
        G4VTransitionRadiation.hh
        G4VXTRenergyLoss.hh
        G4XTRGammaRadModel.hh
        G4XTRRegularRadModel.hh
        G4XTRTransparentRegRadModel.hh
    SOURCES
        G4Cerenkov.cc
        G4ForwardXrayTR.cc
        G4GammaXTRadiator.cc
        G4RegularXTRadiator.cc
        G4Scintillation.cc
        G4ScintillationTrackInformation.cc
        G4StrawTubeXTRadiator.cc
        G4SynchrotronRadiation.cc
        G4SynchrotronRadiationInMat.cc
        G4TransitionRadiation.cc
        G4TransparentRegXTRadiator.cc
        G4VTransitionRadiation.cc
        G4VXTRenergyLoss.cc
        G4XTRGammaRadModel.cc
        G4XTRRegularRadModel.cc
        G4XTRTransparentRegRadModel.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4emutils
        G4geometrymng
        G4globman
        G4hepnumerics
        G4ions
        G4leptons
        G4magneticfield
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

