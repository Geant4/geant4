# - G4emdna-molman module build definition

# Define the Geant4 Module.
geant4_add_module(G4emdna-molman
  PUBLIC_HEADERS
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
    G4MoleculeTableMessenger.hh
  SOURCES
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
    G4MoleculeTableMessenger.cc)

geant4_module_link_libraries(G4emdna-molman
  PUBLIC
    G4emdna-man
    G4globman
    G4partman
  PRIVATE
    G4heprandom
    G4intercoms
    G4track)
