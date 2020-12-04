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
GEANT4_DEFINE_MODULE(NAME G4emdna-molman
    HEADERS
        G4FakeParticleID.hh
        G4MolecularConfiguration.hh
        G4MolecularDissociationChannel.hh
        G4MolecularDissociationTable.hh
        G4MoleculeCounter.hh
        G4MoleculeDefinition.hh
        G4MoleculeFinder.hh
        G4MoleculeHandleManager.hh
        G4Molecule.hh
        G4MoleculeIterator.hh
        G4MoleculeTable.hh
        G4Serialize.hh
        G4VMolecularDissociationDisplacer.hh
        G4VMoleculeCounter.hh
    SOURCES
        G4FakeParticleID.cc
        G4MolecularConfiguration.cc
        G4MolecularDissociationChannel.cc
        G4MolecularDissociationTable.cc
        G4MoleculeCounter.cc
        G4Molecule.cc
        G4MoleculeDefinition.cc
        G4MoleculeHandleManager.cc
        G4MoleculeTable.cc
        G4Serialize.cc
        G4VMolecularDissociationDisplacer.cc
        G4VMoleculeCounter.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4volumes
        G4globman
        G4heprandom
        G4intercoms
        G4materials
        G4partman
        G4procman
        G4track
        G4emdna-man
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

