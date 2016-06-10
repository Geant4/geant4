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
# $Id: sources.cmake 70601 2013-06-03 11:20:53Z gcosmo $
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
include(Geant4MacroDefineModule)


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
    #  !!! ONLY QT4 to be comment with QT5
    include(${QT_USE_FILE})
    
    # Must enable Qt source with a def...
    GEANT4_ADD_COMPILE_DEFINITIONS(SOURCES G4Qt.cc
        COMPILE_DEFINITIONS G4INTY_BUILD_QT)

    # It uses Qt core and gui(?) libs
    list(APPEND G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES
        "${QT_QTGUI_LIBRARY};${QT_QTCORE_LIBRARY}") 
endif()

#
# Add Wt if required
#
if(GEANT4_USE_WT)
    # Add the sources
    list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Wt.hh)
    list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Wt.cc)
        
    # Must have Wt includes...
    include_directories(${Wt_INCLUDE_DIR})

    # Source needs to have a compile definition
    GEANT4_ADD_COMPILE_DEFINITIONS(
        SOURCES G4Wt.cc
        COMPILE_DEFINITIONS G4INTY_BUILD_WT
        )

    # It uses Wt libs
    list(APPEND G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES
        "${Wt_LIBRARY}") 
endif()

#
# Win32 option - always built when we're using MSVC
#
if(MSVC)
    # Add the sources
    list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Win32.hh)
    list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Win32.cc)

    # Add the needed compile definitions to these sources
    GEANT4_ADD_COMPILE_DEFINITIONS(
        SOURCES G4Win32.cc
        COMPILE_DEFINITIONS G4INTY_BUILD_WIN32
    )
endif()

#
# X11 options - we only need this for Xm UI
#
if(UNIX)
  if(GEANT4_USE_XM OR GEANT4_USE_INVENTOR)
    # Add the sources
    list(APPEND G4INTERFACES_COMMON_MODULE_HEADERS G4Xt.hh)
    list(APPEND G4INTERFACES_COMMON_MODULE_SOURCES G4Xt.cc)

    # Must have X11 includes...
    include_directories(${X11_INCLUDE_DIR})

    # Source needs to have a compile definition
    GEANT4_ADD_COMPILE_DEFINITIONS(
        SOURCES G4Xt.cc
        COMPILE_DEFINITIONS G4INTY_BUILD_XT
    )

    # It uses the X11 libs?
    list(APPEND G4INTERFACES_COMMON_MODULE_LINK_LIBRARIES
        "${X11_LIBRARIES}"
    )
  endif()
endif()


#
# Define the Geant4 Module.
#
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

