#------------------------------------------------------------------------------
# sources.cmake
# Module : G4brep
# Package: Geant4.src.G4geometry..G4brep
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.2 2010-09-29 20:47:19 bmorgan Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/Boolean/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4brep
    HEADERS
        G4PlacedSolid.icc
        G4Assembly.hh
        G4Hyperbola.icc
        G4Assembly.icc
        G4PlacedSolid.hh
        G4Axis2Placement3D.hh
        G4PlacementVector.hh
        G4Axis2Placement3D.icc
        G4Plane.hh
        G4BREPSolid.hh
        G4Point3DVector.hh
        G4BREPSolid.icc
        G4PointRat.hh
        G4BREPSolidBox.hh
        G4PointRat.icc
        G4BREPSolidCone.hh
        G4ProjectedSurface.hh
        G4BREPSolidCylinder.hh
        G4Ray.icc
        G4Ray.hh
        G4BREPSolidOpenPCone.hh
        G4Sort.hh
        G4BREPSolidPCone.hh
        G4STEPEntity.hh
        G4BREPSolidPolyhedra.hh
        G4SphericalSurface.hh
        G4BREPSolidSphere.hh
        G4Surface.icc
        G4BREPSolidTorus.hh
        G4Surface.hh
        G4BSplineCurve.hh
        G4SurfaceBoundary.hh
        G4BSplineCurve.icc
        G4SurfaceBoundary.icc
        G4BSplineCurveWithKnots.hh
        G4SurfaceList.hh
        G4BSplineSurface.hh
        G4ThreeMat.hh
        G4BSplineSurface.icc
        G4ToroidalSurface.hh
        G4BSplineSurfaceWithKnots.hh
        G4ToroidalSurface.icc
        G4BezierSurface.hh
        G4UVHit.hh
        G4BezierSurface.icc
        G4BoundingBox3D.hh
        G4BoundingBox3D.icc
        G4CircularCurve.hh
        G4CircularCurve.icc
        G4CompositeCurve.hh
        G4CompositeCurve.icc
        G4Conic.hh
        G4Conic.icc
        G4ConicalSurface.hh
        G4ConicalSurface.icc
        G4ControlPoints.hh
        G4ControlPoints.icc
        G4ConvexHull.hh
        G4Curve.hh
        G4Curve.icc
        G4CurvePoint.hh
        G4CurvePoint.icc
        G4CurveRayIntersection.hh
        G4CurveRayIntersection.icc
        G4CurveVector.hh
        G4CylindricalSurface.hh
        G4CylindricalSurface.icc
        G4Ellipse.hh
        G4Ellipse.icc
        G4FConicalSurface.hh
        G4FConicalSurface.icc
        G4FCylindricalSurface.hh
        G4FCylindricalSurface.icc
        G4FPlane.hh
        G4FPlane.icc
        G4Globals.hh
        G4Hyperbola.hh
        G4KnotVector.hh
        G4KnotVector.icc
        G4Line.hh
        G4Line.icc
        G4OsloMatrix.hh
        G4OsloMatrix.icc
        G4Parabola.hh
        G4Parabola.icc
        G4ProjectedSurface.icc
        G4RectangularTrimmedSurface.hh
        G4SphericalSurface.icc
        G4SurfaceOfLinearExtrusion.hh
        G4SurfaceOfRevolution.hh
    SOURCES
        G4Assembly.cc
        G4Axis2Placement3D.cc
        G4BREPSolid.cc
        G4BREPSolidBox.cc
        G4BREPSolidCone.cc
        G4BREPSolidCylinder.cc
        G4BREPSolidOpenPCone.cc
        G4BREPSolidPCone.cc
        G4BREPSolidPolyhedra.cc
        G4BREPSolidSphere.cc
        G4BREPSolidTorus.cc
        G4BSplineCurve.cc
        G4BSplineCurveWithKnots.cc
        G4BSplineSurface.cc
        G4BSplineSurfaceWithKnots.cc
        G4BezierSurface.cc
        G4BoundingBox3D.cc
        G4CircularCurve.cc
        G4CompositeCurve.cc
        G4Conic.cc
        G4ConicalSurface.cc
        G4ControlPoints.cc
        G4Curve.cc
        G4CurvePoint.cc
        G4CurveRayIntersection.cc
        G4CylindricalSurface.cc
        G4Ellipse.cc
        G4FConicalSurface.cc
        G4FCylindricalSurface.cc
        G4FPlane.cc
        G4Hyperbola.cc
        G4KnotVector.cc
        G4Line.cc
        G4OsloMatrix.cc
        G4Parabola.cc
        G4PlacedSolid.cc
        G4PointRat.cc
        G4ProjectedSurface.cc
        G4Ray.cc
        G4RectangularTrimmedSurface.cc
        G4Sort.cc
        G4SphericalSurface.cc
        G4Surface.cc
        G4SurfaceBoundary.cc
        G4SurfaceList.cc
        G4SurfaceOfLinearExtrusion.cc
        G4SurfaceOfRevolution.cc
        G4ThreeMat.cc
        G4ToroidalSurface.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4geomBoolean
        G4geometrymng
        G4globman
        G4graphics_reps
        G4intercoms
    GLOBAL_DEPENDENCIES
        G4global
        G4graphics_reps
        G4intercoms
    LINK_LIBRARIES
)

# List any source specific properties here

