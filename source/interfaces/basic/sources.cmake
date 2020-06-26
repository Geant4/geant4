#------------------------------------------------------------------------------
# Module : G4UIbasic
# Package: Geant4.src.G4interfaces.G4UIbasic
#------------------------------------------------------------------------------

#
# Module has optional sources
#

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

# Tcsh only on UNIX style systems, but always built here
if(UNIX)
    list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UItcsh.hh)
    list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UItcsh.cc)
endif()

# Win32 terminal only for MSVC builds, but always built here
if(MSVC)
  list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIWin32.hh)
  list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIWin32.cc)

  add_definitions(-DG4UI_BUILD_WIN32_SESSION -DG4INTY_BUILD_WIN32)
endif()

# Qt only if selected.
if(GEANT4_USE_QT)
  list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIQt.hh)
  list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIQt.cc)

  # Add the definitions
  # We have to also add in G4INTY_BUILD_QT 'cause G4QT header needs that...
  add_definitions(-DG4UI_BUILD_QT_SESSION -DG4INTY_BUILD_QT)

  # Add the extra libraries
  list(APPEND G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES Qt5::Gui Qt5::Widgets Qt5::Core)
endif()

#
# Wt (deprecated, left in case it is revived)
#
#if(GEANT4_USE_WT)
#    list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIWt.cc)
#    list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIWt.hh)
#
#    # We have to also add in G4INTY_BUILD_Wt 'cause G4Wt header needs that...
#    GEANT4_ADD_COMPILE_DEFINITIONS(SOURCES G4UIWt.cc G4UIExecutive.cc
#        COMPILE_DEFINITIONS G4UI_BUILD_WT_SESSION;G4INTY_BUILD_WT)
#
#    # Add the extra libraries
#    list(APPEND G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES Wt::Wt Wt::HTTP)
#endif()


# Xm and only on UNIX and if selected
if(UNIX AND GEANT4_USE_XM)
  list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIXm.hh)
  list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIXm.cc)

  # Add the compile definitions - also need INTY versions
  add_definitions(-DG4UI_BUILD_XM_SESSION -DG4INTY_BUILD_XT)

  # Need the X11 and Motif libraries
  list(APPEND G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES
      Motif::Xm
      X11::SM X11::ICE X11::X11 X11::Xext X11::Xmu X11::Xt
  )
endif()

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4UIbasic
  HEADERS
    ${G4INTERFACES_BASIC_MODULE_HEADERS}
  SOURCES
    ${G4INTERFACES_BASIC_MODULE_SOURCES}
  GRANULAR_DEPENDENCIES
    G4UIcommon
    G4globman
    G4intercoms
  GLOBAL_DEPENDENCIES
    G4global
    G4intercoms
  LINK_LIBRARIES
    ${G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES}
    ${timemory_LIBRARIES}
)

# List any source specific properties here
