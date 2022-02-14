# - G4UIbasic module build definition

# Module has optional sources

# List those always built
set(G4INTERFACES_BASIC_MODULE_HEADERS
  G4UIArrayString.hh
  G4UIExecutive.hh
  G4UIcsh.hh
  G4UIterminal.hh
  G4VUIshell.hh)

set(G4INTERFACES_BASIC_MODULE_SOURCES
  G4UIArrayString.cc
  G4UIExecutive.cc
  G4UIcsh.cc
  G4UIterminal.cc
  G4VUIshell.cc)


set(G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES )
set(G4INTERFACES_BASIC_MODULE_COMPILE_DEFINITIONS )

# Tcsh only on UNIX style systems, but always built here
if(UNIX)
  list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UItcsh.hh)
  list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UItcsh.cc)
endif()

# Win32 terminal only for MSVC builds, but always built here
if(MSVC)
  list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIWin32.hh)
  list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIWin32.cc)

  # Definition required by G4UIexecutive
  list(APPEND G4INTERFACES_BASIC_MODULE_COMPILE_DEFINITIONS G4UI_BUILD_WIN32_SESSION)
endif()

# Qt only if selected.
if(GEANT4_USE_QT)
  list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIQt.hh)
  list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIQt.cc)

  # Add the definition required by G4UIexecutive
  list(APPEND G4INTERFACES_BASIC_MODULE_COMPILE_DEFINITIONS G4UI_BUILD_QT_SESSION)

  # Add the extra libraries
  list(APPEND G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES Qt5::Gui Qt5::Widgets Qt5::Core)
endif()

# Xm and only on UNIX and if selected
if(UNIX AND GEANT4_USE_XM)
  list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIXm.hh)
  list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIXm.cc)

  # Add the compile definition needed by G4UIexecutive
  list(APPEND G4INTERFACES_BASIC_MODULE_COMPILE_DEFINITIONS G4UI_BUILD_XM_SESSION)

  # Need the X11 and Motif libraries
  list(APPEND G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES
      Motif::Xm
      X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu X11::Xt
  )
endif()

# Define the Geant4 Module.
geant4_add_module(G4UIbasic
  PUBLIC_HEADERS
    ${G4INTERFACES_BASIC_MODULE_HEADERS}
  SOURCES
    ${G4INTERFACES_BASIC_MODULE_SOURCES})

geant4_module_compile_definitions(G4UIbasic PRIVATE ${G4INTERFACES_BASIC_MODULE_COMPILE_DEFINITIONS})

geant4_module_link_libraries(G4UIbasic
  PUBLIC
    G4UIcommon
    G4intercoms
    G4globman
    ${G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES}
    ${timemory_LIBRARIES})

# List any source specific properties here
if(MSVC)
  # Ensure that the CommCtrl library is loaded in Windows
  geant4_module_link_libraries(G4UIbasic PUBLIC COMCTL32)
endif()
