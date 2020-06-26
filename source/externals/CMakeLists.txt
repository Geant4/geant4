#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4externals
# Package: Geant4.src.G4externals
#
# Intermediate level CMakeLists.txt - just process subdirectories
#
# Created on : 02/06/2011
#
# $Date$
# $Revision$
# $Author$
#
#------------------------------------------------------------------------------

#----------------------------------------------------------------------------
# We only need internal CLHEP if it was requested
#
if(NOT GEANT4_USE_SYSTEM_CLHEP)
    add_subdirectory(clhep)
endif()

#----------------------------------------------------------------------------
# We only need internal expat if it was requested
#
if(NOT GEANT4_USE_SYSTEM_EXPAT)
    add_subdirectory(expat)
endif()

#----------------------------------------------------------------------------
# Geant4 internal ZLIB build optional
#
if(NOT GEANT4_USE_SYSTEM_ZLIB)
    add_subdirectory(zlib)
endif()

#----------------------------------------------------------------------------
# Geant4 internal PTL build optional
#
if(NOT GEANT4_USE_SYSTEM_PTL)
    set(CMAKE_INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}-${${PROJECT_NAME}_VERSION}
        CACHE PATH "CMake install directory" FORCE)
    add_subdirectory(ptl)
    # Needs to be in defined category list for static/shared to be linked correctly.
    set_property(GLOBAL APPEND PROPERTY GEANT4_DEFINED_CATEGORIES G4ptl)
endif()

#----------------------------------------------------------------------------
# Alias the tasking libraries to Geant4 target naming convention
#
if(TARGET ptl-shared)
    add_library(G4ptl INTERFACE)
    add_library(PTL::ptl-shared ALIAS ptl-shared)
    target_link_libraries(G4ptl INTERFACE PTL::ptl-shared)
    install(TARGETS G4ptl
      EXPORT Geant4LibraryDepends
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
    add_library(Geant4::G4ptl ALIAS G4ptl)
endif()

if(TARGET ptl-static)
    add_library(G4ptl-static INTERFACE)
    add_library(PTL::ptl-static ALIAS ptl-static)
    target_link_libraries(G4ptl-static INTERFACE PTL::ptl-static)
    install(TARGETS G4ptl-static
      EXPORT Geant4LibraryDepends
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
    add_library(Geant4::G4ptl-static ALIAS G4ptl-static)
endif()
