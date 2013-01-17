#-----------------------------------------------------------------------
# CMakeLists.txt
# Module : G4global
# Package: Geant4.src.G4global
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id$
#
#-----------------------------------------------------------------------

# Add allocation export symbol for the global category
add_definitions(-DG4GLOB_ALLOC_EXPORT)

add_subdirectory(HEPGeometry)
add_subdirectory(HEPNumerics)
add_subdirectory(HEPRandom)
add_subdirectory(management)

# Geant4 MT libraries - built in isolation, but we ensure they get
# exported into the build tree
# - tpmalloc 
#set(TPMALLOC_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}")
#set(TPMALLOC_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/Geant4")
#set(TPMALLOC_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
#set(TPMALLOC_INSTALL_EXPORT_NAME "Geant4LibraryDepends")
#add_subdirectory(tpmalloc)
#set_property(GLOBAL APPEND
#  PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/source/global/tpmalloc)
#set_property(GLOBAL APPEND
#  PROPERTY GEANT4_EXPORTED_TARGETS tpmallocstub mymalloc)

# - memoryprotection
#set(MEMORYPROTECTION_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}")
#set(MEMORYPROTECTION_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/Geant4")
#set(MEMORYPROTECTION_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
#set(MEMORYPROTECTION_INSTALL_EXPORT_NAME "Geant4LibraryDepends")
#add_subdirectory(memoryprotection)
#set_property(GLOBAL APPEND
#  PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/source/global/memoryprotection)
#set_property(GLOBAL APPEND
#  PROPERTY GEANT4_EXPORTED_TARGETS ProtectedMemory)

# Geant4 global library
if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4global
        COMPONENTS
            HEPGeometry/sources.cmake
            HEPNumerics/sources.cmake
            HEPRandom/sources.cmake
            management/sources.cmake
)
endif()

