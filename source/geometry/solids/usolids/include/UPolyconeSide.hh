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
// UPolyconeSide
//
// Class description:
//
//   Class implmenting a face that represents one conical side
//   of a polycone:
//
//   UPolyconeSide( const UPolyconeSideRZ *prevRZ,
//                   const UPolyconeSideRZ *tail,
//                   const UPolyconeSideRZ *head,
//                   const UPolyconeSideRZ *nextRZ,
//                         double phiStart, double deltaPhi,
//                         bool phiIsOpen, bool isAllBehind=false )
//
//   Values for r1,z1 and r2,z2 should be specified in clockwise
//   order in (r,z).
//
// 19.04.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UPolyconeSide_hh
#define UPolyconeSide_hh

#include "UVCSGface.hh"

class UIntersectingCone;

struct UPolyconeSideRZ
{
  double r, z;  // start of vector
};

class UPolyconeSidePrivateSubclass
{
  public:
    std::pair<UVector3, double> fPhi; // Cached value for phi

    void initialize()
    {
      fPhi.first = UVector3(0, 0, 0);
      fPhi.second = 0.0;
    };
};

class UPolyconeSide : public UVCSGface
{
  public:

    UPolyconeSide(const UPolyconeSideRZ* prevRZ,
                  const UPolyconeSideRZ* tail,
                  const UPolyconeSideRZ* head,
                  const UPolyconeSideRZ* nextRZ,
                  double phiStart, double deltaPhi,
                  bool phiIsOpen, bool isAllBehind = false);
    virtual ~UPolyconeSide();

    UPolyconeSide(const UPolyconeSide& source);
    UPolyconeSide& operator=(const UPolyconeSide& source);

    bool Distance(const UVector3& p, const UVector3& v,
                  bool outgoing, double surfTolerance,
                  double& distance, double& distFromSurface,
                  UVector3& normal, bool& isAllBehind);

    double Safety(const UVector3& p, bool outgoing);

    VUSolid::EnumInside Inside(const UVector3& p, double tolerance,
                               double* bestDistance);

    UVector3 Normal(const UVector3& p,  double* bestDistance);

    double Extent(const UVector3 axis);

    /*
    void CalculateExtent( const EAxisType axis,
                          const UVoxelLimits &voxelLimit,
                          const UAffineTransform &tranform,
                                USolidExtentList &extentList       );
                                */

    UVCSGface* Clone()
    {
      return new UPolyconeSide(*this);
    }

    double SurfaceArea();
    UVector3 GetPointOnFace();

  public: // without description

    UPolyconeSide(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  protected:

    double DistanceAway(const UVector3& p, bool opposite,
                        double& distOutside2, double* rzNorm = 0);

    double DistanceAway(const UVector3& p,
                        double& distOutside2, double* edgeRZnorm);

    bool PointOnCone(const UVector3& hit, double normSign,
                     const UVector3& p,
                     const UVector3& v, UVector3& normal);

    void CopyStuff(const UPolyconeSide& source);

    static void FindLineIntersect(double x1, double y1,
                                  double tx1, double ty1,
                                  double x2, double y2,
                                  double tx2, double ty2,
                                  double& x, double& y);

    double GetPhi(const UVector3& p);

  protected:

    double r[2], z[2]; // r, z parameters, in specified order
    double startPhi,   // Start phi (0 to 2pi), if phiIsOpen
           deltaPhi;   // Delta phi (0 to 2pi), if phiIsOpen
    bool   phiIsOpen; // True if there is a phi slice
    bool   allBehind; // True if the entire solid is "behind" this face

    UIntersectingCone* cone;  // Our intersecting utility class

    double rNorm, zNorm;  // Normal to surface in r,z space
    double rS, zS;        // Unit vector along surface in r,z space
    double length;        // Length of face in r,z space
    double prevRS,
           prevZS;        // Unit vector along previous polyconeSide
    double nextRS,
           nextZS;        // Unit vector along next polyconeSide

    double rNormEdge[2],
           zNormEdge[2];  // Normal to edges

    int ncorners;
    UVector3* corners; // The coordinates of the corners (if phiIsOpen)

  private:
    double tolerance; // Geometrical surface thickness
    double fSurfaceArea;  // Used for surface calculation
};

#endif
