# - G4UIcommon module build definition

# Module has optional sources

# List those always built
set(G4INTERFACES_COMMON_MODULE_HEADERS
  G4InteractorMessenger.hh
  G4VBasicShell.hh
  G4VInteractiveSession.hh
  G4VInteractorManager.hh)

set(G4INTERFACES_COMMON_MODULE_SOURCES
  G4InteractorMessenger.cc
  G4VBasicShell.cc
  G4VInteractiveSession.cc
  G4VInteractorManager.cc)

set(G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES )
set(G4INTERFACES_COMMON_MODULE_COMPILE_DEFINITIONS )

# Add Qt if required
if(GEANT4_USE_QT)
  list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Qt.hh)
  list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Qt.cc)

  # It uses Qt core and gui(?) libs
  list(APPEND G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES
    Qt5::Gui Qt5::Core)
endif()

# Win32 option - always built when we're using MSVC
if(MSVC)
  # Add the sources
  list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Win32.hh)
  list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Win32.cc)
endif()

# X11 options - we only need this for Xm (and Inventor, which activates XM)
if(UNIX AND GEANT4_USE_XM)
  # Add the sources
  list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Xt.hh)
  list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Xt.cc)

  # It uses the X11 libs?
  list(APPEND G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES
    X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu X11::Xt
  )
endif()

if(GEANT4_USE_VTK)
  list(APPEND G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES
       ${VTK_LIBRARIES})
endif()

# Define the Geant4 Module.
geant4_add_module(G4UIcommon
  PUBLIC_HEADERS
    ${G4INTERFACES_COMMON_MODULE_HEADERS}
  SOURCES
    ${G4INTERFACES_COMMON_MODULE_SOURCES})

geant4_module_link_libraries(G4UIcommon
  PUBLIC
    G4intercoms
    G4globman
    # These are PUBLIC until use in client code better understood
    ${G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES})
