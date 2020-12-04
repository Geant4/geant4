#------------------------------------------------------------------------------
# sources.cmake
# Module : G4visQt3D
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
GEANT4_DEFINE_MODULE(NAME G4visQt3D
    HEADERS
        G4Qt3D.hh
        G4Qt3DQEntity.hh
        G4Qt3DSceneHandler.hh
        G4Qt3DUtils.hh
        G4Qt3DViewer.hh
    SOURCES
        G4Qt3D.cc
        G4Qt3DSceneHandler.cc
        G4Qt3DUtils.cc
        G4Qt3DViewer.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hits
        G4intercoms
        G4interfaces
        G4modeling
        G4specsolids
        G4tracking
        G4vis_management
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4geometry
        G4global
        G4graphics_reps
        G4intercoms
        G4interfaces
        G4modeling
        G4tracking
        G4vis_management
    LINK_LIBRARIES
        Qt5::Gui Qt5::Widgets Qt5::3DCore Qt5::3DExtras Qt5::3DRender
)

add_definitions(-DG4INTY_BUILD_QT -DG4UI_BUILD_QT_SESSION -DG4VIS_BUILD_QT3D_DRIVER)

# List any source specific properties here

