#------------------------------------------------------------------------------
# sources.cmake
# Module : G4empii
# Package: Geant4.src.G4processes.G4electromagnetic.G4empii
#
# Sources description for a library.
# Lists the sources and headers of the code explicitly.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 19/11/2010
#
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4empii
    HEADERS
        G4CompositeDataSet.hh
        G4DataSet.hh
        G4hImpactIonisation.hh
        G4hRDEnergyLoss.hh
        G4IDataSet.hh
        G4IInterpolator.hh
        G4LinInterpolator.hh
        G4LogLogInterpolator.hh
        G4PixeCrossSectionHandler.hh
        G4PixeShellDataSet.hh
    SOURCES
        G4CompositeDataSet.cc
        G4DataSet.cc
        G4hImpactIonisation.cc
        G4hRDEnergyLoss.cc
        G4LinInterpolator.cc
        G4LogLogInterpolator.cc
        G4PixeCrossSectionHandler.cc
        G4PixeShellDataSet.cc
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

