# - G4optical module build definition

# Define the Geant4 Module.
geant4_add_module(G4optical
  PUBLIC_HEADERS
    G4OpAbsorption.hh
    G4OpBoundaryProcess.hh
    G4OpMieHG.hh
    G4OpProcessSubType.hh
    G4OpRayleigh.hh
    G4OpWLS.hh
    G4OpWLS2.hh
    G4VWLSTimeGeneratorProfile.hh
    G4WLSTimeGeneratorProfileDelta.hh
    G4WLSTimeGeneratorProfileExponential.hh
  SOURCES
    G4OpAbsorption.cc
    G4OpBoundaryProcess.cc
    G4OpMieHG.cc
    G4OpRayleigh.cc
    G4OpWLS.cc
    G4OpWLS2.cc
    G4VWLSTimeGeneratorProfile.cc
    G4WLSTimeGeneratorProfileDelta.cc
    G4WLSTimeGeneratorProfileExponential.cc)

geant4_module_link_libraries(G4optical
  PUBLIC
    G4bosons
    G4globman
    G4heprandom
    G4materials
    G4procman
  PRIVATE
    G4detector
    G4emutils
    G4navigation
    G4scoring
    G4volumes)
