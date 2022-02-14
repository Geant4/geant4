# - G4solidstate_channeling module build definition

# Define the Geant4 Module.
geant4_add_module(G4solidstate_channeling
  PUBLIC_HEADERS
    G4ChannelingTrackData.hh
    G4ChannelingOptrMultiParticleChangeCrossSection.hh
    G4ChannelingOptrChangeCrossSection.hh
    G4ChannelingMaterialData.hh
    G4Channeling.hh
    G4ChannelingECHARM.hh
  SOURCES
    G4ChannelingTrackData.cc
    G4ChannelingOptrMultiParticleChangeCrossSection.cc
    G4ChannelingOptrChangeCrossSection.cc
    G4ChannelingMaterialData.cc
    G4Channeling.cc
    G4ChannelingECHARM.cc)

geant4_module_link_libraries(G4solidstate_channeling
  PUBLIC
    G4biasing_mgt
    G4geometrymng
    G4globman
    G4materials
    G4procman
    G4track
  PRIVATE
    G4baryons
    G4biasing_gen
    G4emutils
    G4heprandom
    G4partman
    G4volumes)
