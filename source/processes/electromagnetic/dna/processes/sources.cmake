#------------------------------------------------------------------------------
# sources.cmake
# Module : G4emlowenergy
# Package: Geant4.src.G4processes.G4electromagnetic.G4emlowenergy
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 91584 2015-07-27 13:01:48Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/standard/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/lowenergy/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/models/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/molecules/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/molecules/types/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4emdna-processes
    HEADERS
        G4DNAAttachment.hh
        G4DNABrownianTransportation.hh
        G4DNAChargeDecrease.hh
        G4DNAChargeIncrease.hh
        G4DNADissociation.hh
        G4DNAElastic.hh
        G4DNAElectronSolvatation.hh
        G4DNAElectronHoleRecombination.hh
        G4DNAExcitation.hh
        G4DNAIonisation.hh
        G4DNAWaterDissociationDisplacer.hh
        G4DNAMolecularDissociation.hh
        G4DNAPositronium.hh
        G4DNARotExcitation.hh
        G4DNASecondOrderReaction.hh
        G4DNAVibExcitation.hh
    SOURCES
        G4DNAAttachment.cc
        G4DNABrownianTransportation.cc
        G4DNAChargeDecrease.cc
        G4DNAChargeIncrease.cc
        G4DNADissociation.cc
        G4DNAElastic.cc
        G4DNAElectronSolvatation.cc
        G4DNAElectronHoleRecombination.cc
        G4DNAExcitation.cc
        G4DNAIonisation.cc
        G4DNAMolecularDissociation.cc
        G4DNAPositronium.cc
        G4DNARotExcitation.cc
        G4DNAWaterDissociationDisplacer.cc
        G4DNASecondOrderReaction.cc
        G4DNAVibExcitation.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4emlowenergy
        G4emstandard
        G4emutils
        G4geometrymng
        G4navigation
        G4volumes
        G4globman
        G4hepnumerics
        G4intercoms
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4procman
        G4track
        G4emdna-man
        G4emdna-molman
        G4emdna-moltypes
        G4emdna-utils
        G4emdna-models
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

