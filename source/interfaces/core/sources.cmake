# - G4UIcore module build definition

# Define the Geant4 Module.
geant4_add_module(G4UIcore
  PUBLIC_HEADERS
    G4InteractorMessenger.hh
    G4UIArrayString.hh
    G4UIcsh.hh
    G4UIterminal.hh
    G4VBasicShell.hh
    G4VInteractiveSession.hh
    G4VInteractorManager.hh
    G4VUIshell.hh
  SOURCES
    G4InteractorMessenger.cc
    G4UIArrayString.cc
    G4UIcsh.cc
    G4UIterminal.cc
    G4VBasicShell.cc
    G4VInteractiveSession.cc
    G4VInteractorManager.cc
    G4VUIshell.cc)

if(UNIX)
  geant4_module_sources(G4UIcore PUBLIC_HEADERS G4UItcsh.hh SOURCES G4UItcsh.cc)
endif()

geant4_module_link_libraries(G4UIcore PUBLIC G4globman G4intercoms)

