#------------------------------------------------------------------------------
# sources.cmake
# Module : G4visXXX
# Package: Geant4.src.G4visualization.G4visXXX
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
GEANT4_DEFINE_MODULE(NAME G4visXXX
    HEADERS
        G4XXX.hh
        G4XXXFile.hh
        G4XXXFileSceneHandler.hh
        G4XXXFileViewer.hh
        G4XXXSG.hh
        G4XXXSGSceneHandler.hh
        G4XXXSGViewer.hh
        G4XXXSceneHandler.hh
        G4XXXStored.hh
        G4XXXStoredSceneHandler.hh
        G4XXXStoredViewer.hh
        G4XXXViewer.hh
    SOURCES
        G4XXX.cc
        G4XXXFile.cc
        G4XXXFileSceneHandler.cc
        G4XXXFileViewer.cc
        G4XXXSG.cc
        G4XXXSGSceneHandler.cc
        G4XXXSGViewer.cc
        G4XXXSceneHandler.cc
        G4XXXStored.cc
        G4XXXStoredSceneHandler.cc
        G4XXXStoredViewer.cc
        G4XXXViewer.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hits
        G4intercoms
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
        G4modeling
        G4tracking
        G4vis_management
    LINK_LIBRARIES
)

# List any source specific properties here

