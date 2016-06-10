#------------------------------------------------------------------------------
# sources.cmake
# Module : G4geomUSolids
# Package: Geant4.src.G4externals.G4geomUSolids
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

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4geomUSolids
    HEADERS
        UBits.hh
        UBox.hh
        UCons.hh
        UCons.icc
        UEnclosingCylinder.hh
        UExtrudedSolid.hh
        UExtrudedSolid.icc
        UGenericPolycone.hh
        UGenericPolycone.icc
        UGenericTrap.hh
        UGenericTrap.icc
        UIntersectingCone.hh
        UMultiUnion.hh
        UOrb.hh
        UPolycone.hh
        UPolycone.icc
        UPolyconeSide.hh
        UPolyhedra.hh
        UPolyhedra.icc
        UPolyhedraSide.hh
        UPolyPhiFace.hh
        UPolyPhiFace.icc
        UQuadrangularFacet.hh
        UReduciblePolygon.hh
        USphere.hh
        UTessellatedGeometryAlgorithms.hh
        UTessellatedSolid.hh
        UTet.hh
        UTransform3D.hh
        UTrap.hh
        UTrap.icc
        UTrd.hh
        UTrd.icc
        UTriangularFacet.hh
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
        UBits.cc
        UBox.cc
        UCons.cc
        UEnclosingCylinder.cc
        UExtrudedSolid.cc
        UGenericPolycone.cc
        UGenericTrap.cc
        UIntersectingCone.cc
        UMultiUnion.cc
        UOrb.cc
        UPolycone.cc
        UPolyconeSide.cc
        UPolyhedra.cc
        UPolyhedraSide.cc
        UPolyPhiFace.cc
        UQuadrangularFacet.cc
        UReduciblePolygon.cc
        USphere.cc
        UTessellatedGeometryAlgorithms.cc
        UTessellatedSolid.cc
        UTet.cc
        UTransform3D.cc
        UTrap.cc
        UTrd.cc
        UTriangularFacet.cc
        UTubs.cc
        UUtils.cc
        UVCSGfaceted.cc
	UVector2.cc
        UVector3.cc
        UVoxelizer.cc
        VUFacet.cc
        VUSolid.cc
    GRANULAR_DEPENDENCIES
    GLOBAL_DEPENDENCIES
    LINK_LIBRARIES
)

# List any source specific properties here

