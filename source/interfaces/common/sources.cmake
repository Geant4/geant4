#------------------------------------------------------------------------------
# sources.cmake
# Module : G4UIcommon
# Package: Geant4.src.G4interfaces.G4UIcommon
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 18:46:19 bmorgan Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
#include_directories(${CMAKE_SOURCE_DIR}/source/interfaces/basic/include)

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

#
# Add Qt if required
#
if(GEANT4_USE_QT)
    # Add the sources
    list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Qt.hh)
    list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Qt.cc)
    
    # Stuff for Moc etc, plus LINK_LIBRARIES here
endif()

#
# Win32 option
#
if(MSVC AND GEANT4_USE_WIN32TERMINAL)
    # Add the sources
    list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Win32.hh)
    list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Win32.cc)

    # Any extra things here
endif()

#
# X11/Xt options
#
if(UNIX AND GEANT4_USE_X11TERMINAL)
    # Add the sources
    list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Xt.hh)
    list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Xt.cc)

    # Any extra things here
endif()


#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4UIcommon
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

