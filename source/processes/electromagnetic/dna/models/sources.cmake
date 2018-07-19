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
# $Id: sources.cmake 103929 2017-05-03 13:47:48Z gcosmo $
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/molecules/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/molecules/types/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/dna/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4emdna-models
    HEADERS
        G4DNABornAngle.hh
        G4DNABornExcitationModel.hh
        G4DNABornExcitationModel1.hh
        G4DNABornExcitationModel2.hh
        G4DNABornIonisationModel.hh
        G4DNABornIonisationModel1.hh
        G4DNABornIonisationModel2.hh
        G4DNAChampionElasticModel.hh
        G4DNACPA100ElasticModel.hh
        G4DNACPA100ExcitationModel.hh
        G4DNACPA100IonisationModel.hh
        G4DNASmoluchowskiDiffusion.hh
        G4DNASmoluchowskiReactionModel.hh
        G4DNADingfelderChargeDecreaseModel.hh
        G4DNADingfelderChargeIncreaseModel.hh
        G4DNAEmfietzoglouExcitationModel.hh
        G4DNAEmfietzoglouIonisationModel.hh
        G4DNAIonElasticModel.hh
        G4DNAMeltonAttachmentModel.hh
        G4DNAMillerGreenExcitationModel.hh
        G4DNAMolecularReaction.hh
        G4DNAMolecularStepByStepModel.hh
        G4DNAMoleculeEncounterStepper.hh
        G4DNARuddAngle.hh
        G4DNARuddIonisationExtendedModel.hh
        G4DNARuddIonisationModel.hh
        G4DNASancheExcitationModel.hh
        G4DNAOneStepThermalizationModel.hh
        G4DNAOneStepThermalizationModel.hpp
        G4DNAPTBIonisationModel.hh
        G4DNAPTBElasticModel.hh
        G4DNAPTBExcitationModel.hh
        G4DNAPTBAugerModel.hh        
        G4DNAScreenedRutherfordElasticModel.hh
        G4DNATransformElectronModel.hh
        G4DNAUeharaScreenedRutherfordElasticModel.hh
        G4DNAVacuumModel.hh
        G4LEPTSAttachmentModel.hh
        G4LEPTSDissociationModel.hh
        G4LEPTSElasticModel.hh
        G4LEPTSIonisationModel.hh
        G4LEPTSPositroniumModel.hh
        G4LEPTSRotExcitationModel.hh
        G4LEPTSVibExcitationModel.hh
        G4VLEPTSModel.hh
        G4LEPTSDiffXS.hh
        G4LEPTSDistribution.hh
        G4LEPTSElossDistr.hh
        G4LEPTSExcitationModel.hh
        G4VDNAModel.hh
        G4DNAModelInterface.hh
        G4DNADummyModel.hh
    SOURCES
        G4DNABornAngle.cc
        G4DNABornExcitationModel1.cc
        G4DNABornExcitationModel2.cc
        G4DNABornIonisationModel1.cc
        G4DNABornIonisationModel2.cc
        G4DNAChampionElasticModel.cc
        G4DNACPA100ElasticModel.cc
        G4DNACPA100ExcitationModel.cc
        G4DNACPA100IonisationModel.cc
        G4DNASmoluchowskiDiffusion.cc
        G4DNASmoluchowskiReactionModel.cc
        G4DNADingfelderChargeDecreaseModel.cc
        G4DNADingfelderChargeIncreaseModel.cc
        G4DNAEmfietzoglouExcitationModel.cc
        G4DNAEmfietzoglouIonisationModel.cc
        G4DNAIonElasticModel.cc
        G4DNAMeltonAttachmentModel.cc
        G4DNAMillerGreenExcitationModel.cc
        G4DNAMolecularReaction.cc
        G4DNAMolecularStepByStepModel.cc
        G4DNAMoleculeEncounterStepper.cc
        G4DNARuddAngle.cc
        G4DNARuddIonisationExtendedModel.cc
        G4DNARuddIonisationModel.cc
        G4DNASancheExcitationModel.cc
        G4DNAOneStepThermalizationModel.cc
        G4DNAPTBIonisationModel.cc
        G4DNAPTBElasticModel.cc
        G4DNAPTBExcitationModel.cc
        G4DNAPTBAugerModel.cc
        G4DNAScreenedRutherfordElasticModel.cc
        G4DNATransformElectronModel.cc
        G4DNAUeharaScreenedRutherfordElasticModel.cc
        G4DNAVacuumModel.cc
        G4LEPTSElossDistr.cc
        G4LEPTSAttachmentModel.cc
        G4LEPTSDissociationModel.cc
        G4LEPTSElasticModel.cc
        G4LEPTSDistribution.cc
        G4LEPTSIonisationModel.cc
        G4LEPTSPositroniumModel.cc
        G4LEPTSRotExcitationModel.cc
        G4LEPTSVibExcitationModel.cc
        G4VLEPTSModel.cc
        G4LEPTSExcitationModel.cc
        G4LEPTSDiffXS.cc
        G4VDNAModel.cc
        G4DNAModelInterface.cc
        G4DNADummyModel.cc
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
        G4emdna-man
        G4emdna-molman
        G4emdna-moltypes
        G4emdna-utils
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


