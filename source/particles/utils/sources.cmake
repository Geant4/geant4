#------------------------------------------------------------------------------
# Module : G4partutils
# Package: Geant4.src.G4particles.G4partutils
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4partutils
  HEADERS
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
    G4VParticlePropertyReporter.cc
  GRANULAR_DEPENDENCIES
    G4geometrymng
    G4globman
    G4intercoms
    G4materials
    G4partman
  GLOBAL_DEPENDENCIES
    G4geometry
    G4global
    G4intercoms
    G4materials
)

# List any source specific properties here
