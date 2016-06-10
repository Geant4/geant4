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
// UTrap
//
// Class description:
//
// A UTrap is a general trapezoid: The faces perpendicular to the
// z planes are trapezia, and their centres are not necessarily on
// a line parallel to the z axis.
//
// Note that of the 11 parameters described below, only 9 are really
// independent - a check for planarity is made in the calculation of the
// equation for each plane. If the planes are not parallel, a call to
// UException is made.
//
//      pDz    Half-length along the z-axis
//      pTheta  Polar angle of the line joining the centres of the faces
//              at -/+pDz
//      pPhi    Azimuthal angle of the line joing the centre of the face at
//              -pDz to the centre of the face at +pDz
//      pDy1    Half-length along y of the face at -pDz
//      pDx1    Half-length along x of the side at y=-pDy1 of the face at -pDz
//      pDx2    Half-length along x of the side at y=+pDy1 of the face at -pDz
//      pAlp1  Angle with respect to the y axis from the centre of the side
//              at y=-pDy1 to the centre at y=+pDy1 of the face at -pDz
//
//      pDy2    Half-length along y of the face at +pDz
//      pDx3    Half-length along x of the side at y=-pDy2 of the face at +pDz
//      pDx4    Half-length along x of the side at y=+pDy2 of the face at +pDz
//      pAlp2  Angle with respect to the y axis from the centre of the side
//              at y=-pDy2 to the centre at y=+pDy2 of the face at +pDz
//
//   Member Data:
//
//      fDz    Half-length along the z axis
//      fTthetaCphi = std::tan(pTheta)*std::cos(pPhi)
//      fTthetaSphi = std::tan(pTheta)*std::sin(pPhi)
//      These combinations are suitable for creation of the trapezoid corners
//
//      fDy1    Half-length along y of the face at -fDz
//      fDx1    Half-length along x of the side at y=-fDy1 of the face at -fDz
//      fDx2    Half-length along x of the side at y=+fDy1 of the face at -fDz
//      fTalpha1   Tan of Angle with respect to the y axis from the centre of
//                 the side at y=-fDy1 to the centre at y=+fDy1 of the face
//                 at -fDz
//
//      fDy2    Half-length along y of the face at +fDz
//      fDx3    Half-length along x of the side at y=-fDy2 of the face at +fDz
//      fDx4    Half-length along x of the side at y=+fDy2 of the face at +fDz
//      fTalpha2   Tan of Angle with respect to the y axis from the centre of
//                 the side at y=-fDy2 to the centre at y=+fDy2 of the face
//                 at +fDz
//
//      UTrapSidePlane fPlanes[4]  Plane equations of the faces not at +/-fDz
//                                 NOTE: order is important !!!
//
// 12.02.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UTrap_HH
#define UTrap_HH

#include "VUSolid.hh"

struct UTrapSidePlane
{
  double a, b, c, d; // Normal Unit vector (a,b,c) and offset (d)
  // => Ax+By+Cz+D=0
};

class UTrap : public VUSolid
{

  public: // with description

    UTrap(const std::string& pName,
          double pDz,
          double pTheta, double pPhi,
          double pDy1, double pDx1, double pDx2,
          double pAlp1,
          double pDy2, double pDx3, double pDx4,
          double pAlp2);
    //
    // The most general constructor for UTrap which prepares plane
    // equations and corner coordinates from parameters

    UTrap(const std::string& pName,
          const UVector3 pt[8]) ;
    //
    // Prepares plane equations and parameters from corner coordinates

    UTrap(const std::string& pName,
          double pZ,
          double pY,
          double pX, double pLTX);
    //
    // Constructor for Right Angular Wedge from STEP (assumes pLTX<=pX)

    UTrap(const std::string& pName,
          double pDx1,  double pDx2,
          double pDy1,  double pDy2,
          double pDz);
    //
    // Constructor for UTrd

    UTrap(const std::string& pName,
          double pDx, double pDy, double pDz,
          double pAlpha, double pTheta, double pPhi);
    //
    // Constructor for UPara

    UTrap(const std::string& pName);
    //
    // Constructor for "nominal" UTrap whose parameters are to be Set
    // by a UVPVParamaterisation later

    virtual ~UTrap() ;
    //
    // Destructor

    // Accessors

    inline double GetZHalfLength()  const;
    inline double GetYHalfLength1() const;
    inline double GetXHalfLength1() const;
    inline double GetXHalfLength2() const;
    inline double GetTanAlpha1()    const;
    inline double GetYHalfLength2() const;
    inline double GetXHalfLength3() const;
    inline double GetXHalfLength4() const;
    inline double GetTanAlpha2()    const;
    //
    // Returns coordinates of Unit vector along straight
    // line joining centers of -/+fDz planes

    inline UTrapSidePlane GetSidePlane(int n) const;
    inline UVector3 GetSymAxis() const;

    // Modifiers

    void SetAllParameters(double pDz,
                          double pTheta,
                          double pPhi,
                          double pDy1,
                          double pDx1,
                          double pDx2,
                          double pAlp1,
                          double pDy2,
                          double pDx3,
                          double pDx4,
                          double pAlp2);

    void SetPlanes(const UVector3 pt[8]);

    // Methods for solid

    inline double Capacity();
    inline double SurfaceArea();

    VUSolid::EnumInside Inside(const UVector3& p) const;

    UVector3 SurfaceNormal(const UVector3& p) const;

    bool Normal(const UVector3& aPoint, UVector3& aNormal) const;

    double DistanceToIn(const UVector3& p, const UVector3& v, 
                        double aPstep = UUtils::kInfinity) const;

    double SafetyFromOutside(const UVector3& p, bool precise = false) const;

    double DistanceToOut(const UVector3& p,
                         const UVector3&  v,
                         UVector3&       aNormalVector,
                         bool&           aConvex,
                         double aPstep = UUtils::kInfinity) const;

    double SafetyFromInside(const UVector3& p, bool precise = false) const;

    UGeometryType GetEntityType() const;

    UVector3 GetPointOnSurface() const;

    VUSolid* Clone() const;
  
    virtual void Extent(UVector3& aMin, UVector3& aMax) const;

    std::ostream& StreamInfo(std::ostream& os) const;

    // Visualisation functions

  public: // without description

    UTrap(const UTrap& rhs);
    UTrap& operator=(const UTrap& rhs);
    // Copy constructor and assignment operator.

    inline double GetThetaCphi() const;
    inline double GetThetaSphi() const;

  protected:  // with description

    bool MakePlanes();
    bool MakePlane(const UVector3& p1,
                   const UVector3& p2,
                   const UVector3& p3,
                   const UVector3& p4,
                   UTrapSidePlane& plane) ;

  private:

    UVector3 ApproxSurfaceNormal(const UVector3& p) const;
    // Algorithm for SurfaceNormal() following the original
    // specification for points not on the surface

    inline double GetFaceArea(const UVector3& p1,
                              const UVector3& p2,
                              const UVector3& p3,
                              const UVector3& p4);
    //
    // Provided four corners of plane in clockwise fashion,
    // it returns the area of finite face

    UVector3 GetPointOnPlane(UVector3 p0, UVector3 p1,
                             UVector3 p2, UVector3 p3,
                             double& area) const;
    //
    // Returns a random point on the surface of one of the faces

    void GetParametersList(int /*aNumber*/, double* /*aArray*/) const {}


    void ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}

  private:

    double fDz, fTthetaCphi, fTthetaSphi;
    double fDy1, fDx1, fDx2, fTalpha1;
    double fDy2, fDx3, fDx4, fTalpha2;
    UTrapSidePlane fPlanes[4];

    double fCubicVolume;
    double fSurfaceArea;

};

#include "UTrap.icc"

#endif
