#------------------------------------------------------------------------------
# sources.cmake
# Module : G4FR
# Package: Geant4.src.G4visualization.G4FR
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 88190 2015-02-02 17:24:54Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/modeling/include)


#
# Module has optional sources
#
include(Geant4MacroDefineModule)

# List those always built
set(G4VIS_DAWN_MODULE_HEADERS
    G4DAWNFILE.hh
    G4DAWNFILESceneHandler.hh
    G4DAWNFILEViewer.hh
    G4FRConst.hh
    G4FRFeatures.hh
    G4FRSceneFunc.icc
    G4FRofstream.hh
    G4VisFeaturesOfDAWNFILE.hh
    G4VisFeaturesOfFukuiRenderer.hh)

set(G4VIS_DAWN_MODULE_SOURCES
    G4DAWNFILE.cc
    G4DAWNFILESceneHandler.cc
    G4DAWNFILEViewer.cc
    G4VisFeaturesOfDAWNFILE.cc
    G4VisFeaturesOfFukuiRenderer.cc)

#
# DAWN Network driver only built if user selected
#
if(GEANT4_USE_NETWORKDAWN)
    list(APPEND G4VIS_DAWN_MODULE_HEADERS
        G4FRClientServer.hh
        G4FukuiRenderer.hh
        G4FukuiRendererSceneHandler.hh
        G4FukuiRendererViewer.hh)

    list(APPEND G4VIS_DAWN_MODULE_SOURCES
        G4FRClientServer.cc
        G4FukuiRenderer.cc
        G4FukuiRendererSceneHandler.cc
        G4FukuiRendererViewer.cc)

    # To activate the Fukui Renderer Network component, we need an
    # extra compile definition
    GEANT4_ADD_COMPILE_DEFINITIONS(SOURCES ${G4VIS_DAWN_MODULE_SOURCES}
        COMPILE_DEFINITIONS G4VIS_BUILD_DAWN_DRIVER)
endif()



#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4FR
    HEADERS
        ${G4VIS_DAWN_MODULE_HEADERS}
    SOURCES
        ${G4VIS_DAWN_MODULE_SOURCES}
    GRANULAR_DEPENDENCIES
        G4csg
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hits
        G4intercoms
        G4modeling
        G4specsolids
        G4vis_management
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4geometry
        G4global
        G4graphics_reps
        G4intercoms
        G4modeling
        G4vis_management
    LINK_LIBRARIES
)

# List any source specific properties here

