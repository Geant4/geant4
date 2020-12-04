#------------------------------------------------------------------------------
# Module : G4UIcommon
# Package: Geant4.src.G4interfaces.G4UIcommon
#------------------------------------------------------------------------------

#
# Module has optional sources
#

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

# Add Qt if required
if(GEANT4_USE_QT)
  # Add the sources
  list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Qt.hh)
  list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Qt.cc)

  # Must enable Qt source with a def...
  add_definitions(-DG4INTY_BUILD_QT)

  # It uses Qt core and gui(?) libs
  list(APPEND G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES
    Qt5::Gui Qt5::Core)
endif()

#
# Wt (deprecated)
#
#if(GEANT4_USE_WT)
#   # Add the sources
#    list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Wt.hh)
#    list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Wt.cc)
#
#    # Source needs to have a compile definition
#    GEANT4_ADD_COMPILE_DEFINITIONS(
#        SOURCES G4Wt.cc
#        COMPILE_DEFINITIONS G4INTY_BUILD_WT
#        )
#
#    # It uses Wt libs
#    list(APPEND G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES
#        Wt::Wt Wt::HTTP)
#endif()

# Win32 option - always built when we're using MSVC
if(MSVC)
  # Add the sources
  list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Win32.hh)
  list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Win32.cc)

  # Add the needed compile definitions to these sources
  add_definitions(-DG4INTY_BUILD_WIN32)
endif()

# X11 options - we only need this for Xm (and Inventor, which activates XM)
if(UNIX AND GEANT4_USE_XM)
  # Add the sources
  list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Xt.hh)
  list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Xt.cc)

  # Source needs to have a compile definition
  add_definitions(-DG4INTY_BUILD_XT)

  # It uses the X11 libs?
  list(APPEND G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES
    X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu X11::Xt
  )
endif()

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4UIcommon
  HEADERS
    ${G4INTERFACES_COMMON_MODULE_HEADERS}
  SOURCES
    ${G4INTERFACES_COMMON_MODULE_SOURCES}
  GRANULAR_DEPENDENCIES
    G4globman
    G4intercoms
  GLOBAL_DEPENDENCIES
    G4global
    G4intercoms
  LINK_LIBRARIES
    ${G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES}
)

# List any source specific properties here
