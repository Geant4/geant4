# - G4UIbasic module build definition

# Define the Geant4 Module.
geant4_add_module(G4UIbasic
  PUBLIC_HEADERS
    G4UIArrayString.hh
    G4UIExecutive.hh
    G4UIcsh.hh
    G4UIterminal.hh
    G4VUIshell.hh
  SOURCES
    G4UIArrayString.cc
    G4UIExecutive.cc
    G4UIcsh.cc
    G4UIterminal.cc
    G4VUIshell.cc)

geant4_module_link_libraries(G4UIbasic
  PUBLIC
    G4UIcommon
    G4intercoms
    G4globman
    ${timemory_LIBRARIES})

# Tcsh only on UNIX style systems, but always built here
if(UNIX)
  geant4_module_sources(G4UIbasic PUBLIC_HEADERS G4UItcsh.hh SOURCES G4UItcsh.cc)
endif()

# Win32 terminal only for MSVC builds, but always built here
if(MSVC)
  geant4_module_sources(G4UIbasic PUBLIC_HEADERS G4UIWin32.hh SOURCES G4UIWin32.cc)
  geant4_module_compile_definitions(G4UIbasic PRIVATE G4UI_BUILD_WIN32_SESSION)
  geant4_module_link_libraries(G4UIbasic PUBLIC COMCTL32)
endif()

# Qt only if selected.
if(GEANT4_USE_QT)
  geant4_module_sources(G4UIbasic PUBLIC_HEADERS G4UIQt.hh SOURCES G4UIQt.cc)
  geant4_module_compile_definitions(G4UIbasic PRIVATE G4UI_BUILD_QT_SESSION)
  geant4_module_link_libraries(G4UIbasic PUBLIC Qt${QT_VERSION_MAJOR}::Gui Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::Core)
endif()

# Xm and only on UNIX and if selected
if(UNIX AND GEANT4_USE_XM)
  geant4_module_sources(G4UIbasic PUBLIC_HEADERS G4UIXm.hh SOURCES G4UIXm.cc)  
  geant4_module_compile_definitions(G4UIbasic PUBLIC G4UI_BUILD_XM_SESSION)
  geant4_module_link_libraries(G4UIbasic 
    PUBLIC
      Motif::Xm
      X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu X11::Xt)
endif()

