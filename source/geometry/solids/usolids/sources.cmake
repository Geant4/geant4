#------------------------------------------------------------------------------
# sources.cmake
# Module : G4geomUSolids
# Package: Geant4.src.G4geometry.G4geomUSolids
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 66356 2012-12-18 09:02:32Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4geomUSolids
    HEADERS
        G4USolid.hh
        UBits.hh
        UBox.hh
        UCons.hh
        UCons.icc
        UEnclosingCylinder.hh
        UGenericPolycone.hh
        UGenericPolycone.icc
        UIntersectingCone.hh
        UOrb.hh
        UPolycone.hh
        UPolycone.icc
        UPolyconeSide.hh
        UPolyhedra.hh
        UPolyhedra.icc
        UPolyhedraSide.hh
        UPolyPhiFace.hh
        UPolyPhiFace.icc
        UReduciblePolygon.hh
        USphere.hh
        UTet.hh
        UTransform3D.hh
        UTrd.hh
        UTrd.icc
        UTubs.hh
        UTubs.icc
        UTypes.hh
        UUtils.hh
        UVCSGface.hh
        UVCSGfaceted.hh
        UVector2.hh
        UVector2.icc
        UVector3.hh
        UVoxelizer.hh
        VUFacet.hh
        VUSolid.hh
    SOURCES
        G4USolid.cc
        UBits.cc
        UBox.cc
        UCons.cc
        UEnclosingCylinder.cc
        UGenericPolycone.cc
        UIntersectingCone.cc
        UOrb.cc
        UPolycone.cc
        UPolyconeSide.cc
        UPolyhedra.cc
        UPolyhedraSide.cc
        UPolyPhiFace.cc
        UReduciblePolygon.cc
        USphere.cc
        UTet.cc
        UTransform3D.cc
        UTrd.cc
        UTubs.cc
        UUtils.cc
        UVCSGfaceted.cc
	UVector2.cc
        UVector3.cc
        UVoxelizer.cc
        VUFacet.cc
        VUSolid.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4graphics_reps
        G4intercoms
        G4volumes
    GLOBAL_DEPENDENCIES
        G4global
        G4graphics_reps
        G4intercoms
    LINK_LIBRARIES
)

# List any source specific properties here

