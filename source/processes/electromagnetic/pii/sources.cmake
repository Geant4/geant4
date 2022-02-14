# - G4empii module build definition

# Define the Geant4 Module.
geant4_add_module(G4empii
  PUBLIC_HEADERS
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
    G4PixeShellDataSet.cc)

geant4_module_link_libraries(G4empii
  PUBLIC
    G4baryons
    G4cuts
    G4emlowenergy
    G4globman
    G4heprandom
    G4leptons
    G4materials
    G4procman
    G4track
  PRIVATE
    G4bosons
    G4emutils
    G4hepnumerics
    G4partman)