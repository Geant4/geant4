# - G4UIimplementation module build definition

# Define the Geant4 Module.
geant4_add_module(G4UIimplementation PUBLIC_HEADERS G4UIExecutive.hh SOURCES G4UIExecutive.cc)
geant4_module_link_libraries(G4UIimplementation PUBLIC G4UIcore PRIVATE G4intercoms G4globman)

# Win32 terminal only for MSVC builds, but always built here
# ... not in G4UIcore because of the requirement to have G4Win32 interactor manager
if(MSVC)
  geant4_module_sources(G4UIimplementation PUBLIC_HEADERS G4UIWin32.hh G4Win32.hh SOURCES G4UIWin32.cc G4Win32.cc)
  geant4_module_compile_definitions(G4UIimplementation
    PUBLIC G4UI_USE_WIN32
    PRIVATE G4UI_BUILD_WIN32_SESSION)
  geant4_module_link_libraries(G4UIimplementation PUBLIC COMCTL32)
endif()

# Qt only if selected.
if(GEANT4_USE_QT)
  geant4_module_sources(G4UIimplementation PUBLIC_HEADERS G4UIQt.hh G4Qt.hh SOURCES G4UIQt.cc G4Qt.cc)
  geant4_module_compile_definitions(G4UIimplementation
    PUBLIC G4UI_USE_QT
    PRIVATE G4UI_BUILD_QT_SESSION)
  geant4_module_link_libraries(G4UIimplementation PUBLIC G4graphics_reps Qt${QT_VERSION_MAJOR}::Gui Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::Core)
  geant4_set_module_property(G4UIimplementation PROPERTY AUTOMOC ON)

  # Coupling...
  if(GEANT4_USE_VTK)
    geant4_module_compile_definitions(G4UIimplementation PRIVATE G4VIS_USE_VTK_QT)
    geant4_module_link_libraries(G4UIimplementation PRIVATE ${VTK_LIBRARIES})
  endif()
endif()

# Xm and only on UNIX and if selected
if(UNIX AND GEANT4_USE_XM)
  geant4_module_sources(G4UIimplementation PUBLIC_HEADERS G4UIXm.hh G4Xt.hh SOURCES G4UIXm.cc G4Xt.cc)
  geant4_module_compile_definitions(G4UIimplementation 
    PUBLIC G4UI_USE_XM
    PRIVATE G4UI_BUILD_XM_SESSION)
  geant4_module_link_libraries(G4UIimplementation 
    PUBLIC
      Motif::Xm
      X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu X11::Xt)
endif()
