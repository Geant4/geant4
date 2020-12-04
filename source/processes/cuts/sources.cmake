#------------------------------------------------------------------------------
# sources.cmake
# Module : G4cuts
# Package: Geant4.src.G4processes.G4cuts
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
GEANT4_DEFINE_MODULE(NAME G4cuts
    HEADERS
        G4MCCIndexConversionTable.hh
        G4MaterialCutsCouple.hh
        G4PhysicsTableHelper.hh
        G4ProductionCuts.hh
        G4ProductionCutsTable.hh
        G4ProductionCutsTableMessenger.hh
        G4RToEConvForElectron.hh
        G4RToEConvForGamma.hh
        G4RToEConvForPositron.hh
        G4RToEConvForProton.hh
        G4VRangeToEnergyConverter.hh
    SOURCES
        G4MCCIndexConversionTable.cc
        G4MaterialCutsCouple.cc
        G4PhysicsTableHelper.cc
        G4ProductionCuts.cc
        G4ProductionCutsTable.cc
        G4ProductionCutsTableMessenger.cc
        G4RToEConvForElectron.cc
        G4RToEConvForGamma.cc
        G4RToEConvForPositron.cc
        G4RToEConvForProton.cc
        G4VRangeToEnergyConverter.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4intercoms
        G4materials
        G4partman
        G4procman
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
    LINK_LIBRARIES
)

# List any source specific properties here

