//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UPolyPhiFace
//
// Class description:
//
//   Definition of a face that bounds a polycone or polyhedra when
//   it has a phi opening:
//
//   UPolyPhiFace( const UReduciblePolygon *rz,
//                        double phi,
//                        double deltaPhi,
//                        double phiOther )
//
//   Specifically: a face that lies on a plane that passes through
//   the z axis. It has boundaries that are straight lines of arbitrary
//   length and direction, but with corners aways on the same side of
//   the z axis.
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UPolyPhiFace_hh
#define UPolyPhiFace_hh

#include "UVCSGface.hh"
#include "UVector2.hh"

class UReduciblePolygon;

struct UPolyPhiFaceVertex
{
  double x, y, r, z;   // position
  double rNorm,
         zNorm;       // r/z normal
  UVector3 norm3D;  // 3D normal

  // Needed for Triangulation Algorithm
  //
  bool ear;
  UPolyPhiFaceVertex* next, *prev;
};

struct UPolyPhiFaceEdge
{
  UPolyPhiFaceEdge(): v0(0), v1(0), tr(.0), tz(0.), length(0.) {}
  UPolyPhiFaceVertex*  v0, *v1; // Corners
  double tr, tz,                // Unit vector along edge
         length;                // Length of edge
  UVector3 norm3D;           // 3D edge normal vector
};

class UPolyPhiFace : public UVCSGface
{

  public: // with description

    UPolyPhiFace(const UReduciblePolygon* rz,
                 double phi, double deltaPhi, double phiOther);
    // Constructor.
    // Points r,z should be supplied in clockwise order in r,z.
    // For example:
    //                [1]---------[2]        ^ R
    //                 |           |          |
    //                 |           |          +--> z
    //                [0]---------[3]

    virtual ~UPolyPhiFace();
    // Destructor. Removes edges and corners.

    UPolyPhiFace(const UPolyPhiFace& source);
    UPolyPhiFace& operator=(const UPolyPhiFace& source);
    // Copy constructor and assgnment operator.

    bool Distance(const UVector3& p, const UVector3& v,
                  bool outgoing, double surfTolerance,
                  double& distance, double& distFromSurface,
                  UVector3& normal, bool& allBehind);

    double Safety(const UVector3& p, bool outgoing);

    VUSolid::EnumInside Inside(const UVector3& p, double tolerance,
                               double* bestDistance);

    UVector3 Normal(const UVector3& p, double* bestDistance);

    double Extent(const UVector3 axis);

    /*
    void CalculateExtent( const EAxisType axis,
                          const UVoxelLimits &voxelLimit,
                          const UAffineTransform &tranform,
                                USolidExtentList &extentList );
                                */

    inline UVCSGface* Clone();
    // Allocates on the heap a clone of this face.

    double SurfaceArea();
    double SurfaceTriangle(UVector3 p1, UVector3 p2,
                           UVector3 p3, UVector3* p4);
    UVector3 GetPointOnFace();
    // Auxiliary methods for determination of points on surface.

  public: // without description

    UPolyPhiFace(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

    void Diagnose(VUSolid* solid);
    // Throw an exception if something is found inconsistent with
    // the solid. For debugging purposes only

  protected:

    bool InsideEdgesExact(double r, double z, double normSign,
                          const UVector3& p, const UVector3& v);
    // Decide if the point in r,z is inside the edges of our face,
    // **but** do so consistently with other faces.

    bool InsideEdges(double r, double z);
    bool InsideEdges(double r, double z, double* distRZ2,
                     UPolyPhiFaceVertex** base3Dnorm = 0,
                     UVector3** head3Dnorm = 0);
    // Decide if the point in r,z is inside the edges of our face.

    inline double ExactZOrder(double z,
                              double qx, double qy, double qz,
                              const UVector3& v,
                              double normSign,
                              const UPolyPhiFaceVertex* vert) const;
    // Decide precisely whether a trajectory passes to the left, right,
    // or exactly passes through the z position of a vertex point in face.

    void CopyStuff(const UPolyPhiFace& source);

  protected:

    // Functions used for Triangulation in Case of generic Polygone.
    // The triangulation is used for GetPointOnFace()

    double Area2(UVector2 a, UVector2 b, UVector2 c);
    // Calculation of 2*Area of Triangle with Sign

    bool Left(UVector2 a, UVector2 b, UVector2 c);
    bool LeftOn(UVector2 a, UVector2 b, UVector2 c);
    bool Collinear(UVector2 a, UVector2 b, UVector2 c);
    // Boolean functions for sign of Surface

    bool IntersectProp(UVector2 a, UVector2 b,
                       UVector2 c, UVector2 d);
    // Boolean function for finding proper intersection of two
    // line segments (a,b) and (c,d).

    bool Between(UVector2 a, UVector2 b, UVector2 c);
    // Boolean function for determining if point c is between a and b
    // where the three points (a,b,c) are on the same line.

    bool Intersect(UVector2 a, UVector2 b,
                   UVector2 c, UVector2 d);
    // Boolean function for finding proper intersection or not
    // of two line segments (a,b) and (c,d).

    bool Diagonalie(UPolyPhiFaceVertex* a, UPolyPhiFaceVertex* b);
    // Boolean Diagonalie help to determine if diagonal s
    // of segment (a,b) is convex or reflex.

    bool InCone(UPolyPhiFaceVertex* a, UPolyPhiFaceVertex* b);
    // Boolean function for determining if b is inside the cone (a0,a,a1)
    // where a is the center of the cone.

    bool Diagonal(UPolyPhiFaceVertex* a, UPolyPhiFaceVertex* b);
    // Boolean function for determining if Diagonal is possible
    // inside Polycone or PolyHedra.

    void EarInit();
    // Initialisation for Triangulisation by ear tips.
    // For details see "Computational Geometry in C" by Joseph O'Rourke.

    void Triangulate();
    // Triangularisation by ear tips for Polycone or Polyhedra.
    // For details see "Computational Geometry in C" by Joseph O'Rourke.
    // NOTE: a copy of the shape is made and this copy is reordered in
    //       order to have a list of triangles. This list is used by the
    //       method GetPointOnFace().

  protected:

    int     numEdges;           // Number of edges
    UPolyPhiFaceEdge*   edges;     // The edges of the face
    UPolyPhiFaceVertex* corners;   // And the corners
    UVector3    normal;       // Normal Unit vector
    UVector3    radial;       // Unit vector along radial direction
    UVector3    surface;       // Point on surface
    UVector3    surface_point; // Auxiliary point on surface used for
    // method GetPointOnFace()
    double   rMin, rMax, // Extent in r
             zMin, zMax; // Extent in z
    bool      allBehind; // True if the polycone/polyhedra
    // is behind the place of this face
    double   fTolerance;// Surface thickness
    double   fSurfaceArea; // Surface Area of PolyPhiFace
    UPolyPhiFaceVertex* triangles; // Auxiliary pointer to 'corners' used for
    // triangulation. Copy structure, changing
    // the structure of 'corners' (ear removal)
};

#include "UPolyPhiFace.icc"

#endif
