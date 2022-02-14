# - G4hadronic_qmd module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_qmd
  PUBLIC_HEADERS
    G4QMDCollision.hh
    G4QMDGroundStateNucleus.hh
    G4QMDMeanField.hh
    G4QMDNucleus.hh
    G4QMDParameters.hh
    G4QMDParticipant.hh
    G4QMDReaction.hh
    G4QMDSystem.hh
  SOURCES
    G4QMDCollision.cc
    G4QMDGroundStateNucleus.cc
    G4QMDMeanField.cc
    G4QMDNucleus.cc
    G4QMDParameters.cc
    G4QMDParticipant.cc
    G4QMDReaction.cc
    G4QMDSystem.cc)

geant4_module_link_libraries(G4hadronic_qmd
  PUBLIC
    G4globman
    G4had_im_r_matrix
    G4hadronic_deex_evaporation
    G4hadronic_deex_handler
    G4hadronic_mgt
    G4hepgeometry
    G4partman
  PRIVATE
    G4baryons
    G4hadronic_util
    G4hadronic_xsect
    G4heprandom
    G4materials)
