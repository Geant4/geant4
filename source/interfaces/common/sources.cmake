# - G4UIcommon module build definition

# Define the Geant4 Module.
geant4_add_module(G4UIcommon
  PUBLIC_HEADERS
    G4InteractorMessenger.hh
    G4VBasicShell.hh
    G4VInteractiveSession.hh
    G4VInteractorManager.hh
  SOURCES
    G4InteractorMessenger.cc
    G4VBasicShell.cc
    G4VInteractiveSession.cc
    G4VInteractorManager.cc)

geant4_module_link_libraries(G4UIcommon PUBLIC G4intercoms G4globman)

# Add Qt if required
if(GEANT4_USE_QT)
  geant4_module_sources(G4UIcommon PUBLIC_HEADERS G4Qt.hh SOURCES G4Qt.cc)
  geant4_module_link_libraries(G4UIcommon PUBLIC Qt${QT_VERSION_MAJOR}::Gui Qt${QT_VERSION_MAJOR}::Core Qt${QT_VERSION_MAJOR}::Widgets)
endif()

# Win32 option - always built when we're using MSVC
if(MSVC)
  geant4_module_sources(G4UIcommon PUBLIC_HEADERS G4Win32.hh SOURCES G4Win32.cc)
endif()

# X11 options - we only need this for Xm (and Inventor, which activates XM)
if(UNIX AND GEANT4_USE_XM)
  geant4_module_sources(G4UIcommon PUBLIC_HEADERS G4Xt.hh SOURCES G4Xt.cc)
  geant4_module_link_libraries(G4UIcommon PUBLIC X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu X11::Xt)
endif()

if(GEANT4_USE_VTK)
  geant4_module_compile_definitions(G4UIcommon PRIVATE G4VIS_USE_VTK_QT)
  geant4_module_link_libraries(G4UIcommon PUBLIC ${VTK_LIBRARIES})
endif()
