# - G4had_string_frag module build definition

# Define the Geant4 Module.
geant4_add_module(G4had_string_frag
  PUBLIC_HEADERS
    G4ExcitedStringDecay.hh
    G4FragmentingString.hh
    G4HadronBuilder.hh
    G4LundStringFragmentation.hh
    G4QGSMFragmentation.hh
    G4VLongitudinalStringDecay.hh
  SOURCES
    G4ExcitedStringDecay.cc
    G4FragmentingString.cc
    G4HadronBuilder.cc
    G4LundStringFragmentation.cc
    G4QGSMFragmentation.cc
    G4VLongitudinalStringDecay.cc)

geant4_module_link_libraries(G4had_string_frag
  PUBLIC
    G4globman
    G4had_string_man
    G4hadronic_mgt
    G4hadronic_util
    G4hepgeometry
    G4partman
  PRIVATE
    G4heprandom
    G4procman
    G4shortlived
    G4track)
