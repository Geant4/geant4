# - G4analysismng module build definition

if(GEANT4_USE_FREETYPE)
  set(G4analysismng_LINK_LIBRARIES Freetype::Freetype)
endif()

# Define the Geant4 Module.
geant4_add_module(G4analysismng
  PUBLIC_HEADERS
    G4AnalysisVerbose.hh
    G4AnalysisManagerState.hh
    G4AnalysisMessenger.hh
    G4AnalysisUtilities.hh
    G4BaseAnalysisManager.hh
    G4BaseFileManager.hh
    G4BaseNtupleManager.hh
    G4BaseRNtupleManager.hh
    G4BinScheme.hh
    G4Fcn.hh
    G4NtupleBookingManager.hh
    G4NtupleBookingManager.icc
    G4NtupleMessenger.hh
    G4HnInformation.hh
    G4HnManager.hh
    G4HnMessenger.hh
    G4TFileInformation.hh
    G4TFileManager.hh
    G4TFileManager.icc
    G4THnManager.hh
    G4THnManager.icc
    G4THnMessenger.hh
    G4THnMessenger.icc
    G4THnToolsManager.hh
    G4THnToolsManager.icc
    G4TNtupleDescription.hh
    G4TNtupleManager.hh
    G4TNtupleManager.icc
    G4TRNtupleDescription.hh
    G4TRNtupleManager.hh
    G4TRNtupleManager.icc
    G4VAnalysisManager.hh
    G4VAnalysisManager.icc
    G4VAnalysisReader.hh
    G4VAnalysisReader.icc
    G4VFileManager.hh
    G4VNtupleFileManager.hh
    G4VNtupleManager.hh
    G4VRFileManager.hh
    G4VRNtupleManager.hh
    G4VTBaseHnManager.hh
    G4VTFileManager.hh
    G4VTHnFileManager.hh
    G4VTHnRFileManager.hh
  SOURCES
    G4AnalysisVerbose.cc
    G4AnalysisManagerState.cc
    G4AnalysisMessenger.cc
    G4AnalysisUtilities.cc
    G4BaseAnalysisManager.cc
    G4BaseFileManager.cc
    G4BaseNtupleManager.cc
    G4BaseRNtupleManager.cc
    G4BinScheme.cc
    G4Fcn.cc
    G4NtupleBookingManager.cc
    G4NtupleMessenger.cc
    G4HnInformation.cc
    G4HnManager.cc
    G4HnMessenger.cc
    G4VAnalysisManager.cc
    G4VAnalysisReader.cc
    G4VFileManager.cc
    G4VNtupleFileManager.cc)

# NB Freetype may be private
geant4_module_link_libraries(G4analysismng PUBLIC G4globman G4intercoms G4tools ${G4analysismng_LINK_LIBRARIES})
