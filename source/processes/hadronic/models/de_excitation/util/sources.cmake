# - G4hadronic_deex_util module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_deex_util
  PUBLIC_HEADERS
    G4AlphaCoulombBarrier.hh
    G4CameronGilbertPairingCorrections.hh
    G4CameronGilbertShellCorrections.hh
    G4CameronShellPlusPairingCorrections.hh
    G4CameronTruranHilfPairingCorrections.hh
    G4CameronTruranHilfShellCorrections.hh
    G4ChatterjeeCrossSection.hh
    G4ConstantLevelDensityParameter.hh
    G4CookPairingCorrections.hh
    G4CookShellCorrections.hh
    G4CoulombBarrier.hh
    G4DeuteronCoulombBarrier.hh
    G4He3CoulombBarrier.hh
    G4KalbachCrossSection.hh
    G4NeutronCoulombBarrier.hh
    G4PairingCorrection.hh
    G4ProtonCoulombBarrier.hh
    G4ShellCorrection.hh
    G4TritonCoulombBarrier.hh
    G4VCoulombBarrier.hh
    G4VLevelDensityParameter.hh
  SOURCES
    G4AlphaCoulombBarrier.cc
    G4CameronGilbertPairingCorrections.cc
    G4CameronGilbertShellCorrections.cc
    G4CameronShellPlusPairingCorrections.cc
    G4CameronTruranHilfPairingCorrections.cc
    G4CameronTruranHilfShellCorrections.cc
    G4ChatterjeeCrossSection.cc
    G4ConstantLevelDensityParameter.cc
    G4CookPairingCorrections.cc
    G4CookShellCorrections.cc
    G4CoulombBarrier.cc
    G4DeuteronCoulombBarrier.cc
    G4He3CoulombBarrier.cc
    G4KalbachCrossSection.cc
    G4NeutronCoulombBarrier.cc
    G4PairingCorrection.cc
    G4ProtonCoulombBarrier.cc
    G4ShellCorrection.cc
    G4TritonCoulombBarrier.cc
    G4VCoulombBarrier.cc)

geant4_module_link_libraries(G4hadronic_deex_util
  PUBLIC
    G4globman
  PRIVATE
    G4hadronic_util)


