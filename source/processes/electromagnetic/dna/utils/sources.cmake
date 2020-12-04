#------------------------------------------------------------------------------
# sources.cmake
# Module : G4emlowenergy
# Package: Geant4.src.G4processes.G4electromagnetic.G4emlowenergy
#
# Sources description for a library.
# Lists the sources and headers of the code explicitly.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
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
	G4ErrorFunction.hh
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
	G4ErrorFunction.cc
        G4MoleculeGun.cc
        G4MoleculeGunMessenger.cc
        G4ReactionTableMessenger.cc
        G4VDNAReactionModel.cc
        G4VUserChemistryList.cc
        # physchemIO
        G4VPhysChemIO.cc
        G4PhysChemIO.cc
    GRANULAR_DEPENDENCIES
        G4analysismng
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

