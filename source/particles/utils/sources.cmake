# - G4partutils module build definition

# Define the Geant4 Module.
geant4_add_module(G4partutils
  PUBLIC_HEADERS
    G4HtmlPPReporter.hh
    G4IsotopeMagneticMomentTable.hh
    G4SimplePPReporter.hh
    G4TextPPReporter.hh
    G4TextPPRetriever.hh
    G4VParticlePropertyReporter.hh
    G4VParticlePropertyRetriever.hh
  SOURCES
    G4HtmlPPReporter.cc
    G4IsotopeMagneticMomentTable.cc
    G4SimplePPReporter.cc
    G4TextPPReporter.cc
    G4TextPPRetriever.cc
    G4VParticlePropertyReporter.cc)

geant4_module_link_libraries(G4partutils PUBLIC G4partman G4globman)
