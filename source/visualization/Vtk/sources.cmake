#------------------------------------------------------------------------------
# sources.cmake
# Module : G4visVtk
# Package: Geant4.src.G4visualization.G4visQt3D
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#

set(G4VIS_MODULE_VTK_HEADERS
    G4Vtk.hh
    G4VtkMessenger.hh
    G4VtkSceneHandler.hh
    G4VtkViewer.hh
    vtkTensorGlyphColor.h)

set(G4VIS_MODULE_VTK_SOURCES
    G4Vtk.cc
    G4VtkMessenger.cc
    G4VtkSceneHandler.cc
    G4VtkViewer.cc
    vtkTensorGlyphColor.cxx)

set(G4VIS_MODULE_VTK_DEFINITIONS
    -DG4VIS_USE_VTK)

if(GEANT4_USE_QT)
    list(APPEND G4VIS_MODULE_VTK_HEADERS
          G4VtkQt.hh
          G4VtkQtSceneHandler.hh
          G4VtkQtViewer.hh)

    list(APPEND G4VIS_MODULE_VTK_SOURCES
         G4VtkQt.cc
         G4VtkQtSceneHandler.cc
         G4VtkQtViewer.cc)

    list(APPEND G4VIS_MODULE_VTK_DEFINITIONS -DG4VIS_USE_VTK_QT)
endif()

add_definitions(${G4VIS_MODULE_VTK_DEFINITIONS})

geant4_add_module(G4visVtk
        PUBLIC_HEADERS ${G4VIS_MODULE_VTK_HEADERS}
        SOURCES ${G4VIS_MODULE_VTK_SOURCES})

geant4_module_compile_definitions(G4visVtk PRIVATE ${G4VIS_MODULE_VTK_DEFINITIONS})
geant4_module_compile_definitions(G4UIcommon PRIVATE ${G4VIS_MODULE_VTK_DEFINITIONS})

geant4_module_link_libraries(G4visVtk
    PUBLIC
        G4modeling
        G4globman
        G4hepgeometry
        G4vis_management
        ${VTK_LIBRARIES}
    PRIVATE
        G4graphics_reps
        G4UIbasic
        G4UIcommon)


