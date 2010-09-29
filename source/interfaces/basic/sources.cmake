#------------------------------------------------------------------------------
# sources.cmake
# Module : G4UIbasic
# Package: Geant4.src.G4interfaces.G4UIbasic
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 18:46:03 bmorgan Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/interfaces/common/include)

#
# Module has optional sources
#
# List those always built
set(G4INTERFACES_BASIC_MODULE_HEADERS 
    G4UIArrayString.hh
    G4UIExecutive.hh
    G4UIExecutive.icc
    G4UIcsh.hh
    G4UIterminal.hh
    G4VUIshell.hh)

set(G4INTERFACES_BASIC_MODULE_SOURCES
    G4UIArrayString.cc
    G4UIcsh.cc
    G4UIterminal.cc
    G4VUIshell.cc)


set(G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES )


#
# Tcsh only on UNIX style systems
#
if(UNIX)
    list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UItcsh.hh)
    list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UItcsh.cc)
endif()


#
# Win32 terminal only on Win32 and if selected.
#
if(WIN32 AND GEANT4_USE_WIN32TERMINAL)
    list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIWin32.hh)
    list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIWin32.cc)
endif()


#
# Qt only if selected.
#
if(GEANT4_USE_QT)
    list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIQt.hh)
    list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIQt.cc)

    # Any other moc things, extra libraries here
endif()


#
# Xm and Xaw only on UNIX and if selected
#
if(UNIX)
    if(GEANT4_USE_XM)
        list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIXm.hh)
        list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIXm.cc)

        # Any other includes, libraries here
    endif()

    if(GEANT4_USE_XAW)
        list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIXaw.hh)
        list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIXaw.cc)

        # Any other includes, libraries here
    endif()
endif()





#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4UIbasic
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
)

# List any source specific properties here

