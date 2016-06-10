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
// UGenericTrap
//
// Class description:
//
// UGenericTrap is a solid which represents an arbitrary trapezoid with
// up to 8 vertices standing on two parallel planes perpendicular to Z axis.
//
// Parameters in the constructor:
// - name               - solid name
// - halfZ              - the solid half length in Z
// - vertices           - the (x,y) coordinates of vertices:
//                        o first four points: vertices[i], i<4
//                          are the vertices sitting on the -halfZ plane;
//                        o last four points: vertices[i], i>=4
//                          are the vertices sitting on the +halfZ plane.
//
//   The order of defining the vertices of the solid is the following:
//      - point 0 is connected with points 1,3,4
//      - point 1 is connected with points 0,2,5
//      - point 2 is connected with points 1,3,6
//      - point 3 is connected with points 0,2,7
//      - point 4 is connected with points 0,5,7
//      - point 5 is connected with points 1,4,6
//      - point 6 is connected with points 2,5,7
//      - point 7 is connected with points 3,4,6
// Points can be identical in order to create shapes with less than
// 8 vertices.
//
// 21.10.13 Tatiana Nikitina, CERN; Ivana Hrivnacova, IPN Orsay
//          Adapted from Root Arb8 implementation
// --------------------------------------------------------------------

#ifndef USOLIDS_UGenericTrap_HH
#define USOLIDS_UGenericTrap_HH

#ifndef USOLIDS_VUSolid
#include "VUSolid.hh"
#endif

#ifndef USOLIDS_UUtils
#include "UUtils.hh"
#endif

#include <vector>

#include "UVector2.hh"

class VUFacet;
class UTessellatedSolid;
class UBox;

class UGenericTrap : public VUSolid
{
  public:  // with description

    UGenericTrap(const std::string& name, double halfZ,
                 const std::vector<UVector2>& vertices);
    // Constructor

    ~UGenericTrap();
    // Destructor

    // Accessors

    inline double    GetZHalfLength() const;
    inline void      SetZHalfLength(double);
    inline int       GetNofVertices() const;
    inline UVector2  GetVertex(int index) const;
    inline const std::vector<UVector2>& GetVertices() const;
    inline double    GetTwistAngle(int index) const;
    inline bool      IsTwisted() const;
    inline int       GetVisSubdivisions() const;
    inline void      SetVisSubdivisions(int subdiv);

    // Solid methods

    EnumInside Inside(const UVector3& aPoint) const;
    bool Normal(const UVector3& aPoint, UVector3& aNormal) const;
    double  SafetyFromInside(const UVector3& aPoint,
                             bool aAccurate = false) const;
    double  SafetyFromOutside(const UVector3& aPoint,
                              bool aAccurate = false) const;
    double  DistanceToIn(const UVector3& aPoint,
                         const UVector3& aDirection,
                         double aPstep = UUtils::kInfinity) const;

    double DistanceToOut(const UVector3& aPoint,
                         const UVector3& aDirection,
                         UVector3&       aNormalVector,
                         bool&           aConvex,
                         double aPstep = UUtils::kInfinity) const;
    void Extent(UVector3& aMin, UVector3& aMax) const;
    double Capacity() ;
    double SurfaceArea() ;
    VUSolid* Clone() const ;

    inline UGeometryType GetEntityType() const { return "GenericTrap"; }
    inline void ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}
    inline void GetParametersList(int /*aNumber*/, double* /*aArray*/) const {}

    UVector3 GetPointOnSurface() const;

    std::ostream& StreamInfo(std::ostream& os) const;

  public:

    UGenericTrap();
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    UGenericTrap(const UGenericTrap& rhs);
    UGenericTrap& operator=(const UGenericTrap& rhs); 
      // Copy constructor and assignment operator.

    void Initialise(const std::vector<UVector2>&  vertices);
    inline UVector3 GetMinimumBBox() const;
    inline UVector3 GetMaximumBBox() const;

  private:

    // Internal methods

    inline void SetTwistAngle(int index, double twist);
    bool  ComputeIsTwisted() ;
    bool  CheckOrder(const std::vector<UVector2>& vertices) const;
    bool  IsSegCrossing(const UVector2& a, const UVector2& b,
                        const UVector2& c, const UVector2& d) const;
    bool  IsSegCrossingZ(const UVector2& a, const UVector2& b,
                         const UVector2& c, const UVector2& d) const;
    bool IsSameLineSegment(const UVector2& p,
                           const UVector2& l1, const UVector2& l2) const;
    bool IsSameLine(const UVector2& p,
                    const UVector2& l1, const UVector2& l2) const;

    void ReorderVertices(std::vector<UVector3>& vertices) const;
    void ComputeBBox();

    VUFacet* MakeDownFacet(const std::vector<UVector3>& fromVertices,
                           int ind1, int ind2, int ind3) const;
    VUFacet* MakeUpFacet(const std::vector<UVector3>& fromVertices,
                         int ind1, int ind2, int ind3) const;
    VUFacet* MakeSideFacet(const UVector3& downVertex0,
                           const UVector3& downVertex1,
                           const UVector3& upVertex1,
                           const UVector3& upVertex0) const;
    UTessellatedSolid* CreateTessellatedSolid() const;

    EnumInside InsidePolygone(const UVector3& p,
                              const UVector2* poly)const;
    double DistToPlane(const UVector3& p,
                       const UVector3& v, const int ipl) const ;
    double DistToTriangle(const UVector3& p,
                          const UVector3& v, const int ipl) const;
    UVector3 NormalToPlane(const UVector3& p,
                           const int ipl) const;
    double SafetyToFace(const UVector3& p, const int iseg) const;
    double GetFaceSurfaceArea(const UVector3& p0,
                              const UVector3& p1,
                              const UVector3& p2,
                              const UVector3& p3) const;
  private:

    // static data members

    static const int       fgkNofVertices;
    static const double    fgkTolerance;

    // data members

    double                 fDz;
    std::vector<UVector2>  fVertices;
    bool                   fIsTwisted;
    double                 fTwist[4];
    UTessellatedSolid*     fTessellatedSolid;
    UVector3               fMinBBoxVector;
    UVector3               fMaxBBoxVector;
    int                    fVisSubdivisions;
    UBox*                  fBoundBox;

    enum ESide {kUndefined, kXY0, kXY1, kXY2, kXY3, kMZ, kPZ};
    // Codes for faces (kXY[num]=num of lateral face,kMZ= minus z face etc)

    double                 fSurfaceArea;
    double                 fCubicVolume;
    // Surface and Volume
};

#include "UGenericTrap.icc"

#endif
