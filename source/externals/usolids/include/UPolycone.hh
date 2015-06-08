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
// UPolycone
//
// Class description:
//
//   Class implementing a CSG-like type "PCON".
//
//   UPolycone( const std::string& name, 
//              double phiStart,     // initial phi starting angle
//              double phiTotal,     // total phi angle
//              int numZPlanes,      // number of z planes
//              const double zPlane[],  // position of z planes
//              const double rInner[],  // tangent distance to inner surface
//              const double rOuter[])  // tangent distance to outer surface
//
//   Alternative constructor, but limited to increasing-only Z sections:
//
//   UPolycone( const std::string& name, 
//              double phiStart,   // initial phi starting angle
//              double phiTotal,   // total phi angle
//              int    numRZ,      // number corners in r,z space
//              const double r[],  // r coordinate of these corners
//              const double z[])  // z coordinate of these corners
//
// 19.04.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UPolycone_hh
#define UPolycone_hh

#include "VUSolid.hh"

#include "UPolyconeSide.hh"
#include "UVCSGfaceted.hh"
#include "UVoxelizer.hh"

#include "UCons.hh"
#include "UTubs.hh"
#include "UBox.hh"

class UEnclosingCylinder;
class UReduciblePolygon;
class UPolyconeHistorical
{
  public:
    UPolyconeHistorical();
    ~UPolyconeHistorical();
    UPolyconeHistorical(const UPolyconeHistorical& source);
    UPolyconeHistorical& operator=(const UPolyconeHistorical& right);

    double fStartAngle;
    double fOpeningAngle;
    int  fNumZPlanes;
    std::vector<double> fZValues;
    std::vector<double> Rmin;
    std::vector<double> Rmax;
};

class UPolycone : public VUSolid
{

  public: // with description

    void Init(
      double phiStart,        // initial phi starting angle
      double phiTotal,        // total phi angle
      int numZPlanes,         // number of z planes
      const double zPlane[],  // position of z planes
      const double rInner[],  // tangent distance to inner surface
      const double rOuter[]);

    UPolycone(const std::string& name) : VUSolid(name)
    {
    }

    UPolycone(const std::string& name,
              double phiStart,        // initial phi starting angle
              double phiTotal,        // total phi angle
              int numZPlanes,         // number of z planes
              const double zPlane[],  // position of z planes
              const double rInner[],  // tangent distance to inner surface
              const double rOuter[]); // tangent distance to outer surface


    UPolycone(const std::string& name,
              double phiStart,    // initial phi starting angle
              double phiTotal,    // total phi angle
              int   numRZ,        // number corners in r,z space
              const double r[],   // r coordinate of these corners
              const double z[]);  // z coordinate of these corners


    virtual ~UPolycone();

    void Reset();

//  inline void SetOriginalParameters(UPolyconeHistorical* pars);

//  inline void SetOriginalParameters();

    std::ostream& StreamInfo(std::ostream& os) const;

    VUSolid::EnumInside Inside(const UVector3& p) const;

    double DistanceToIn(const UVector3& p, const UVector3& v, double aPstep = UUtils::kInfinity) const;

    double  SafetyFromInside(const UVector3& aPoint,
                             bool aAccurate = false) const;
    double  SafetyFromOutside(const UVector3& aPoint,
                              bool aAccurate = false) const;

    double DistanceToOut(const UVector3& aPoint,
                         const UVector3& aDirection,
                         UVector3&       aNormalVector,
                         bool&           aConvex,
                         double aPstep = UUtils::kInfinity) const;

    bool Normal(const UVector3& aPoint, UVector3& aNormal) const;
    //  virtual void Extent ( EAxisType aAxis, double &aMin, double &aMax ) const;
    void Extent(UVector3& aMin, UVector3& aMax) const;
    double Capacity();
    double SurfaceArea();
    UGeometryType GetEntityType() const;

    void ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}

    // Visualisation
    void GetParametersList(int /*aNumber*/, double* /*aArray*/) const {}
    VUSolid* Clone() const;

    UPolycone(const UPolycone& source);
    UPolycone& operator=(const UPolycone& source);
    // Copy constructor and assignment operator.
    void CopyStuff(const UPolycone& source);
    UVector3 GetPointOnSurface() const;

// Methods for random point generation

    UVector3 GetPointOnCone(double fRmin1, double fRmax1,
                            double fRmin2, double fRmax2,
                            double zOne,   double zTwo,
                            double& totArea) const;

    UVector3 GetPointOnTubs(double fRMin, double fRMax,
                            double zOne,  double zTwo,
                            double& totArea) const;

    UVector3 GetPointOnCut(double fRMin1, double fRMax1,
                           double fRMin2, double fRMax2,
                           double zOne,   double zTwo,
                           double& totArea) const;

    UVector3 GetPointOnRing(double fRMin, double fRMax,
                            double fRMin2, double fRMax2,
                            double zOne) const;

    inline double GetStartPhi() const
    {
      return startPhi;
    }

    inline double GetEndPhi() const
    {
      return endPhi;
    }

    inline bool IsOpen() const
    {
      return phiIsOpen;
    }

    inline bool IsGeneric() const
    {
      return false;
    }

    inline int GetNumRZCorner() const
    {
      return numCorner;
    }

    inline UPolyconeSideRZ GetCorner(int index) const
    {
      return corners[index];
    }

    inline UPolyconeHistorical* GetOriginalParameters() const
    {
      return fOriginalParameters;
    }

    inline void SetOriginalParameters(UPolyconeHistorical* pars)
    {
      if (!pars)
        // UException("UPolycone3::SetOriginalParameters()", "GeomSolids0002",
        //            FatalException, "NULL pointer to parameters!");
        *fOriginalParameters = *pars;
    }

  protected:  // without description

//  int fNumSides;
    bool SetOriginalParameters(UReduciblePolygon* rz);
    // Here are our parameters

    double startPhi;    // Starting phi value (0 < phiStart < 2pi)
    double endPhi;      // end phi value (0 < endPhi-phiStart < 2pi)
    bool   phiIsOpen;  // true if there is a phi segment
    int  numCorner;   // number RZ points
    UPolyconeSideRZ* corners; // corner r,z points
    UPolyconeHistorical* fOriginalParameters; // original input parameters
    double       fCubicVolume;   // Cubic Volume
    double       fSurfaceArea;   // Surface Area
    mutable UBox fBox;  // Bounding box for Polycone

    inline void SetOriginalParameters()
    {
      int numPlanes = (int)numCorner / 2;

      fOriginalParameters = new UPolyconeHistorical;

      fOriginalParameters->fZValues.resize(numPlanes);
      fOriginalParameters->Rmin.resize(numPlanes);
      fOriginalParameters->Rmax.resize(numPlanes);

      for (int j = 0; j < numPlanes; j++)
      {
        fOriginalParameters->fZValues[j] = corners[numPlanes + j].z;
        fOriginalParameters->Rmax[j] = corners[numPlanes + j].r;
        fOriginalParameters->Rmin[j] = corners[numPlanes - 1 - j].r;
      }

      fOriginalParameters->fStartAngle = startPhi;
      fOriginalParameters->fOpeningAngle = endPhi - startPhi;
      fOriginalParameters->fNumZPlanes = numPlanes;
    }

    UEnclosingCylinder* enclosingCylinder;

    struct UPolyconeSection
    {
      VUSolid* solid;// true if all points in section are concave in regards to whole polycone, will be determined
      double shift;
      bool tubular;
//    double left, right;
      bool convex; // TURE if all points in section are concave in regards to whole polycone, will be determined, currently not implemented
    };

    std::vector<double> fZs; // z coordinates of given sections
    std::vector<UPolyconeSection> fSections;
    int fMaxSection;

    inline VUSolid::EnumInside InsideSection(int index, const UVector3& p) const;

    inline double SafetyFromInsideSection(int index, const double rho,
                                          const UVector3& p) const
    {
      const UPolyconeSection& section = fSections[index];
      UVector3 ps(p.x(), p.y(), p.z() - section.shift);
      double res=0;
      if (section.tubular)
      {
        UTubs* tubs = (UTubs*) section.solid;
        res = tubs->SafetyFromInsideR(ps,rho, true);
      }
      else
      {
        UCons* cons = (UCons*) section.solid;
        res = cons->SafetyFromInsideR(ps,rho, true);
      }
      return res;
    }

    // Auxiliary method used in SafetyFromInside for finding safety
    // from section in R and Phi
    //
    inline double SafetyFromOutsideSection(int index, const double rho,
                                           const UVector3& p) const
    {
      const UPolyconeSection& section = fSections[index];
      UVector3 ps(p.x(), p.y(), p.z());
      double res=0;
      if (section.tubular)
      {
        UTubs* tubs = (UTubs*) section.solid;
        res = tubs->SafetyFromOutsideR(ps,rho, true);
      }
      else
      {
        UCons* cons = (UCons*) section.solid;
        res = cons->SafetyFromOutsideR(ps,rho, true);
      }
      return res;
    }

    // Auxiliary method used in SafetyFromOutside for finding safety
    // from section
    //
    inline double SafetyFromOutsideSection(int index, const UVector3& p) const
    {
      const UPolyconeSection& section = fSections[index];
      UVector3 ps(p.x(), p.y(),p.z() - section.shift);
      double res=0;
     
      res = section.solid->SafetyFromOutside(ps, true);
      return res;
    }

    bool NormalSection(int index, const UVector3& p, UVector3& n) const
    {
      const UPolyconeSection& section = fSections[index];
      UVector3 ps(p.x(), p.y(), p.z() - section.shift);
      bool res = section.solid->Normal(ps, n);
      return res;
    }

    inline int GetSection(double z) const
    {
      int section = UVoxelizer::BinarySearch(fZs, z);
      if (section < 0) section = 0;
      else if (section > fMaxSection) section = fMaxSection;
      return section;
    }
};

#endif
