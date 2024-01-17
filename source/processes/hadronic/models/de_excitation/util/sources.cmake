# - G4hadronic_deex_util module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_deex_util
  PUBLIC_HEADERS
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
    G4KalbachCrossSection.hh
    G4PairingCorrection.hh
    G4ShellCorrection.hh
    G4VCoulombBarrier.hh
    G4VLevelDensityParameter.hh
  SOURCES
    G4CameronGilbertPairingCorrections.cc
    G4CameronGilbertShellCorrections.cc
    G4CameronShellPlusPairingCorrections.cc
    G4CameronTruranHilfPairingCorrections.cc
    G4CameronTruranHilfShellCorrections.cc
    G4ChatterjeeCrossSection.cc
    G4CookPairingCorrections.cc
    G4CookShellCorrections.cc
    G4CoulombBarrier.cc
    G4KalbachCrossSection.cc
    G4PairingCorrection.cc
    G4ShellCorrection.cc
    G4VCoulombBarrier.cc)

geant4_module_link_libraries(G4hadronic_deex_util
  PUBLIC
    G4globman
  PRIVATE
    G4hadronic_util)


