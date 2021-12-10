# - G4csv module build definition

# Define the Geant4 Module.
geant4_add_module(G4csv
  PUBLIC_HEADERS
    G4CsvAnalysisManager.hh
    G4CsvAnalysisManager.icc
    G4CsvAnalysisReader.hh
    G4CsvAnalysisReader.icc
    G4CsvFileManager.hh
    G4CsvHnFileManager.hh
    G4CsvHnFileManager.icc
    G4CsvHnRFileManager.hh
    G4CsvHnRFileManager.icc
    G4CsvNtupleFileManager.hh
    G4CsvNtupleManager.hh
    G4CsvRFileManager.hh
    G4CsvRNtupleManager.hh
    g4csv_defs.hh
  SOURCES
    G4CsvAnalysisManager.cc
    G4CsvAnalysisReader.cc
    G4CsvFileManager.cc
    G4CsvNtupleFileManager.cc
    G4CsvNtupleManager.cc
    G4CsvRFileManager.cc
    G4CsvRNtupleManager.cc)

geant4_module_link_libraries(G4csv PUBLIC G4analysismng G4hntools G4globman G4tools)
