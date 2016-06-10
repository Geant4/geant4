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
// UPolyhedra
//
// Class description:
//
//   Class implementing a CSG-like type "PGON":
//
//   UPolyhedra( const std::string& name,
//                double phiStart,         - initial phi starting angle
//                double phiTotal,         - total phi angle
//                int numSide,             - number sides
//                int numZPlanes,          - number of z planes
//                const double zPlane[],   - position of z planes
//                const double rInner[],   - tangent distance to inner surface
//                const double rOuter[]  ) - tangent distance to outer surface
//
//   UPolyhedra( const std::string& name,
//                double phiStart,    - initial phi starting angle
//                double phiTotal,    - total phi angle
//                int    numSide,     - number sides
//                int    numRZ,       - number corners in r,z space
//                const double r[],   - r coordinate of these corners
//                const double z[] )  - z coordinate of these corners
//
// 19.09.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UPolyhedra_hh
#define UPolyhedra_hh

#include "UVCSGfaceted.hh"
#include "UPolyhedraSide.hh"

class UEnclosingCylinder;
class UReduciblePolygon;
class UPolyhedraHistorical
{
  public:

    UPolyhedraHistorical();
    ~UPolyhedraHistorical();
    UPolyhedraHistorical(const UPolyhedraHistorical& source);
    UPolyhedraHistorical& operator=(const UPolyhedraHistorical& right);

    double fStartAngle;
    double fOpeningAngle;
    int fNumSide;
    int fNumZPlanes;
    std::vector<double> fZValues;
    std::vector<double> Rmin;
    std::vector<double> Rmax;
};

class UPolyhedra : public UVCSGfaceted
{
  protected:

    inline UPolyhedra(const std::string& name) : UVCSGfaceted(name) {}

  public:  // with description

    void Init(
      double phiStart,    // initial phi starting angle
      double phiTotal,    // total phi angle
      int numSide,        // number sides
      int numZPlanes,     // number of z planes
      const double zPlane[],    // position of z planes
      const double rInner[],    // tangent distance to inner surface
      const double rOuter[]);   // tangent distance to outer surface

    UPolyhedra(const std::string& name,
               double phiStart,    // initial phi starting angle
               double phiTotal,    // total phi angle
               int numSide,        // number sides
               int numZPlanes,     // number of z planes
               const double zPlane[],    // position of z planes
               const double rInner[],    // tangent distance to inner surface
               const double rOuter[]);   // tangent distance to outer surface

    UPolyhedra(const std::string& name,
               double phiStart,    // initial phi starting angle
               double phiTotal,    // total phi angle
               int    numSide,     // number sides
               int    numRZ,       // number corners in r,z space
               const double r[],         // r coordinate of these corners
               const double z[]);        // z coordinate of these corners

    virtual ~UPolyhedra();

    // Methods for solid

    void GetParametersList(int /*aNumber*/, double* /*aArray*/) const {}

    void ComputeBBox(UBBox* /*aBox*/, bool /*aStore*/)
    {
      // Computes bounding box.
      std::cout << "ComputeBBox - Not implemented" << std::endl;
    }

    VUSolid::EnumInside Inside(const UVector3& p) const;

    // double DistanceToInDelete( const UVector3 &p,
    //                      const UVector3 &v ) const;

    double SafetyFromOutside(const UVector3& aPoint, bool aAccurate = false) const;

    UGeometryType  GetEntityType() const;

    VUSolid* Clone() const;

    UVector3 GetPointOnSurface() const;

    std::ostream& StreamInfo(std::ostream& os) const;

    bool Reset();

    // Accessors

    inline int GetNumSide()     const;
    inline double GetStartPhi() const;
    inline double GetEndPhi()   const;
    inline bool IsOpen()        const;
    inline bool IsGeneric()     const;
    inline int GetNumRZCorner() const;
    inline UPolyhedraSideRZ GetCorner(const int index) const;

    inline UPolyhedraHistorical* GetOriginalParameters();
    // Returns internal scaled parameters.
    inline void SetOriginalParameters(UPolyhedraHistorical& pars);
    // Sets internal parameters. Parameters 'Rmin' and 'Rmax' in input must
    // be scaled first by a factor computed as 'cos(0.5*phiTotal/theNumSide)',
    // if not already scaled.

  public:  // without description

    double DistanceToIn(const UVector3& p,
                        const UVector3& v, double aPstep = UUtils::kInfinity) const;

    UPolyhedra(const UPolyhedra& source);
    UPolyhedra& operator=(const UPolyhedra& source);
    // Copy constructor and assignment operator.

    void Extent(UVector3& aMin, UVector3& aMax) const;

  protected:  // without description

    inline void SetOriginalParameters();
    // Sets internal parameters for the generic constructor.

    void Create(double phiStart,            // initial phi starting angle
                double phiTotal,           // total phi angle
                int    numSide,            // number sides
                UReduciblePolygon* rz);     // rz coordinates
    // Generates the shape and is called by each constructor, after the
    // conversion of the arguments

    void CopyStuff(const UPolyhedra& source);
    void DeleteStuff();

    // Methods for generation of random points on surface

    UVector3 GetPointOnPlane(UVector3 p0, UVector3 p1,
                             UVector3 p2, UVector3 p3) const;
    UVector3 GetPointOnTriangle(UVector3 p0, UVector3 p1,
                                UVector3 p2) const;
    UVector3 GetPointOnSurfaceCorners() const;

    
  protected:  // without description

    int   fNumSides;      // Number of sides
    double fStartPhi;    // Starting phi value (0 < phiStart < 2pi)
    double fEndPhi;      // end phi value (0 < endPhi-phiStart < 2pi)
    bool   fPhiIsOpen;   // true if there is a phi segment
    bool   fGenericPgon; // true if created through the 2nd generic constructor
    int   fNumCorner;    // number RZ points
    UPolyhedraSideRZ* fCorners;  // our corners
    UPolyhedraHistorical fOriginalParameters;  // original input parameters
    UEnclosingCylinder* fEnclosingCylinder;

};

#include "UPolyhedra.icc"

#endif
