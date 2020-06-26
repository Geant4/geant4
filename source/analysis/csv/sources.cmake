#------------------------------------------------------------------------------
# Module : G4csv
# Package: Geant4.src.G4analysis.G4csv
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4csv
  HEADERS
    G4CsvAnalysisManager.hh
    G4CsvAnalysisManager.icc
    G4CsvAnalysisReader.hh
    G4CsvAnalysisReader.icc
    G4CsvFileManager.hh
    G4CsvNtupleManager.hh
    G4CsvRFileManager.hh
    G4CsvRNtupleManager.hh
    g4csv_defs.hh
    g4csv.hh
  SOURCES
    G4CsvAnalysisManager.cc
    G4CsvAnalysisReader.cc
    G4CsvFileManager.cc
    G4CsvNtupleManager.cc
    G4CsvRFileManager.cc
    G4CsvRNtupleManager.cc
  GRANULAR_DEPENDENCIES
    G4globman
    G4intercoms
    G4analysismng
    G4hntools
  GLOBAL_DEPENDENCIES
    G4global
    G4intercoms
)

# List any source specific properties here
