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
# $Id: sources.cmake 75000 2013-10-25 10:58:17Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/interfaces/common/include)
include_directories(${CMAKE_SOURCE_DIR}/source/interfaces/GAG/include)

#
# Module has optional sources
#
include(Geant4MacroDefineModule)

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


#
# Tcsh only on UNIX style systems, but always built here
#
if(UNIX)
    list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UItcsh.hh)
    list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UItcsh.cc)
endif()


#
# Win32 terminal only for MSVC builds, but always built here
#
if(MSVC)
    list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIWin32.hh)
    list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIWin32.cc)

    GEANT4_ADD_COMPILE_DEFINITIONS(
        SOURCES G4UIWin32.cc G4UIExecutive.cc
        COMPILE_DEFINITIONS G4UI_BUILD_WIN32_SESSION;G4INTY_BUILD_WIN32
    )
endif()


#
# Qt only if selected.
#
if(GEANT4_USE_QT)
    list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIQt.hh)
    list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIQt.cc)

    # Now need to add Qt in to build moc wrappers
    include(${QT_USE_FILE})

    # Add the moc sources - must use an absolute path to the files being
    # wrapped
    QT4_WRAP_CPP(G4INTERFACES_MOC_SOURCES 
        ${CMAKE_SOURCE_DIR}/source/interfaces/basic/include/G4UIQt.hh 
        OPTIONS -DG4UI_BUILD_QT_SESSION)

    list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES ${G4INTERFACES_MOC_SOURCES})

    # Add the definitions
    # We have to also add in G4INTY_BUILD_QT 'cause G4QT header needs that...
    GEANT4_ADD_COMPILE_DEFINITIONS(SOURCES G4UIQt.cc G4UIExecutive.cc
        COMPILE_DEFINITIONS G4UI_BUILD_QT_SESSION;G4INTY_BUILD_QT)

    # and for the moc file...
    set_source_files_properties(${G4INTERFACES_MOC_SOURCES}
        PROPERTIES COMPILE_DEFINITIONS G4UI_BUILD_QT_SESSION)

    # Add the extra libraries - seem to need to quote variables to get
    # both linked in.
    list(APPEND G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES 
        "${QT_QTGUI_LIBRARY};${QT_QTCORE_LIBRARY}") 
endif()


#
# Wt only if selected.
#
if(GEANT4_USE_WT)
    list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIWt.hh)
    list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIWt.cc)

    # Must have Wt includes...
    include_directories(${Wt_INCLUDE_DIR})

    # Add the definitions
    # We have to also add in G4INTY_BUILD_Wt 'cause G4Wt header needs that...
    GEANT4_ADD_COMPILE_DEFINITIONS(SOURCES G4UIWt.cc G4UIExecutive.cc
        COMPILE_DEFINITIONS G4UI_BUILD_WT_SESSION;G4INTY_BUILD_WT)

    # Add the extra libraries - seem to need to quote variables to get
    # both linked in.
    list(APPEND G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES 
        "${Wt_LIBRARY}") 
endif()


#
# Xm and only on UNIX and if selected
#
if(UNIX)
    if(GEANT4_USE_XM)
        list(APPEND G4INTERFACES_BASIC_MODULE_HEADERS G4UIXm.hh)
        list(APPEND G4INTERFACES_BASIC_MODULE_SOURCES G4UIXm.cc)

        # Need the X11 and XM headers
        include_directories(${X11_INCLUDE_DIR})
        include_directories(${MOTIF_INCLUDE_DIR})

        # Add the compile definitions - also need INTY versions
        GEANT4_ADD_COMPILE_DEFINITIONS(
            SOURCES G4UIXm.cc G4UIExecutive.cc
            COMPILE_DEFINITIONS G4UI_BUILD_XM_SESSION;G4INTY_BUILD_XT
        )

        # Need the X11 and Motif libraries
        list(APPEND G4INTERFACES_BASIC_MODULE_LINK_LIBRARIES
            "${MOTIF_LIBRARIES};${X11_LIBRARIES}"
        )
    endif()
endif()


#
# Define the Geant4 Module.
#
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

