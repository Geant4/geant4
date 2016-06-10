#------------------------------------------------------------------------------
# sources.cmake
# Module : G4graphics_reps
# Package: Geant4.src.G4graphics_reps
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 66376 2012-12-18 09:42:59Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4graphics_reps
    HEADERS
        G4AttCheck.hh
        G4AttDef.hh
        G4AttDefStore.hh
        G4AttDefT.hh
        G4AttHolder.hh
        G4AttUtils.hh
        G4AttValue.hh
        G4Circle.hh
        G4Color.hh
        G4Colour.hh
        G4ConversionFatalError.hh
        G4ConversionUtils.hh
        G4CreatorFactoryT.hh
        G4DimensionedDouble.hh
        G4DimensionedThreeVector.hh
        G4DimensionedType.hh
        G4PlacedPolyhedron.hh
        G4Point3DList.hh
        G4Polyhedron.hh
        G4PolyhedronArbitrary.hh
        G4Polyline.hh
        G4Polymarker.hh
        G4Polymarker.icc
        G4Scale.hh
        G4Scale.icc
        G4SmartFilter.hh
        G4Square.hh
        G4Square.icc
        G4Text.hh
        G4Text.icc
        G4TypeKey.hh
        G4TypeKeyT.hh
        G4VFilter.hh
        G4VGraphicsScene.hh
        G4VMarker.hh
        G4VMarker.icc
        G4VVisManager.hh
        G4VisAttributes.hh
        G4VisAttributes.icc
        G4VisExtent.hh
        G4Visible.hh
        G4Visible.icc
        HepPolyhedron.h
        HepPolyhedronProcessor.h
    SOURCES
        BooleanProcessor.src
        G4AttCheck.cc
        G4AttDef.cc
        G4AttDefStore.cc
        G4AttHolder.cc
        G4AttUtils.cc
        G4Circle.cc
        G4Colour.cc
        G4DimensionedTypeUtils.cc
        G4PlacedPolyhedron.cc
        G4Point3DList.cc
        G4Polyhedron.cc
        G4PolyhedronArbitrary.cc
        G4Polyline.cc
        G4Polymarker.cc
        G4Scale.cc
        G4Square.cc
        G4Text.cc
        G4VGraphicsScene.cc
        G4VMarker.cc
        G4VVisManager.cc
        G4VisAttributes.cc
        G4VisExtent.cc
        G4Visible.cc
        HepPolyhedron.cc
        HepPolyhedronProcessor.src
    GRANULAR_DEPENDENCIES
        G4globman
        G4intercoms
    GLOBAL_DEPENDENCIES
        G4global
        G4intercoms
    LINK_LIBRARIES
)

# List any source specific properties here

