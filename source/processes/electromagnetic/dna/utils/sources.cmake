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
# $Id: sources.cmake 103042 2017-03-10 11:50:07Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/molecules/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/molecules/types/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

include_directories(${CMAKE_SOURCE_DIR}/source/analysis/management/include) # for physchemIO

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4emdna-utils
    HEADERS
        G4DNAChemistryManager.hh
        G4DNACPA100LogLogInterpolation.hh
        G4DNACPA100WaterExcitationStructure.hh
        G4DNACPA100WaterIonisationStructure.hh
        G4DNACrossSectionDataSet.hh
        G4DNADamage.hh
        G4DNAGenericIonsManager.hh
        G4DNAIons.hh
        G4DNAMolecularMaterial.hh
        G4DNAMolecularReactionTable.hh
        G4DNAEmfietzoglouWaterExcitationStructure.hh
        G4DNAEmfietzoglouWaterIonisationStructure.hh
        G4DNAPTBIonisationStructure.hh
        G4DNARevertProbability.hh
        G4DNAWaterExcitationStructure.hh
        G4DNAWaterIonisationStructure.hh
        G4MoleculeGun.hh
        G4MoleculeGunMessenger.hh
        G4ReactionTableMessenger.hh
        G4VDNAReactionModel.hh
        G4VUserChemistryList.hh
        # physchemIO
        G4VPhysChemIO.hh
        G4PhysChemIO.hh
    SOURCES
        G4DNAChemistryManager.cc
        G4DNACPA100LogLogInterpolation.cc
        G4DNACPA100WaterExcitationStructure.cc
        G4DNACPA100WaterIonisationStructure.cc
        G4DNACrossSectionDataSet.cc
        G4DNADamage.cc
        G4DNAGenericIonsManager.cc
        G4DNAIons.cc
        G4DNAMolecularMaterial.cc
        G4DNAMolecularReactionTable.cc
        G4DNAEmfietzoglouWaterExcitationStructure.cc
        G4DNAEmfietzoglouWaterIonisationStructure.cc
        G4DNAPTBIonisationStructure.cc
        G4DNAWaterExcitationStructure.cc
        G4DNAWaterIonisationStructure.cc
        G4MoleculeGun.cc
        G4MoleculeGunMessenger.cc
        G4ReactionTableMessenger.cc
        G4VDNAReactionModel.cc
        G4VUserChemistryList.cc
        # physchemIO
        G4VPhysChemIO.cc
        G4PhysChemIO.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4emlowenergy
        G4emstandard
        G4emutils
        G4geometrymng
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
        # G4emdna-man
        # G4emdna-molman
        # G4emdna-moltypes
    GLOBAL_DEPENDENCIES
        G4analysis
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

