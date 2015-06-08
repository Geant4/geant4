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
// UGenericPolycone
//
// Class description:
//
//   Implementing a CSG-like type "PCON" volume with possibility of
//   specifying also 'decreasing' Z sections:
//
//   UGenericPolycone( const std::string& name,
//                     double phiStart,  // initial phi starting angle
//                     double phiTotal,  // total phi angle
//                     int    numRZ,      // number corners in r,z space
//                     const double r[],  // r coordinate of these corners
//                     const double z[])  // z coordinate of these corners
//
// 19.10.13 Tatiana Nikitina
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UGenericPolycone_hh
#define UGenericPolycone_hh

#include "UVCSGfaceted.hh"
#include "UPolyconeSide.hh"

class UEnclosingCylinder;
class UReduciblePolygon;
class UVCSGface;

class UGenericPolycone: public UVCSGfaceted
{

  public:  // with description

    UGenericPolycone(const std::string& name,
                     double phiStart,     // initial phi starting angle
                     double phiTotal,     // total phi angle
                     int numZPlanes,     // number of z planes
                     const double zPlane[],     // position of z planes
                     const double rInner[],     // tangent distance to inner surface
                     const double rOuter[]);   // tangent distance to outer surface

    UGenericPolycone(const std::string& name,
                     double phiStart,    // initial phi starting angle
                     double phiTotal,    // total phi angle
                     int   numRZ,       // number corners in r,z space
                     const double r[],        // r coordinate of these corners
                     const double z[]);        // z coordinate of these corners

    virtual ~UGenericPolycone();

    // Methods for solid

    VUSolid::EnumInside Inside(const UVector3& p) const;
    double DistanceToIn(const UVector3& p, const UVector3& v, double aPstep = UUtils::kInfinity) const;
    // double SafetyFromOutside( const UVector3 &p, bool aAccurate=false) const;

    UVector3 GetPointOnSurface() const;

    /*
    void ComputeDimensions(      UVPVParameterisation* p,
                            const int n,
                            const UVPhysicalVolume* pRep );
                            */

    UGeometryType GetEntityType() const;

    VUSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;
    bool Reset();

    // Accessors

    inline double GetStartPhi() const;
    inline double GetEndPhi()   const;
    inline bool IsOpen()         const;
    inline int  GetNumRZCorner() const;
    inline UPolyconeSideRZ GetCorner(int index) const;


  public:  // without description

    //UPolycone(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

    UGenericPolycone(const UGenericPolycone& source);
    UGenericPolycone& operator=(const UGenericPolycone& source);
    // Copy constructor and assignment operator.

  protected: // without description

    // Generic initializer, called by all constructors


    void Create(double phiStart,         // initial phi starting angle
                double phiTotal,       // total phi angle
                UReduciblePolygon* rz);  // r/z coordinate of these corners

    void CopyStuff(const UGenericPolycone& source);

    // Methods for random point generation



    void GetParametersList(int /*aNumber*/, double* /*aArray*/) const {}

    void ComputeBBox(UBBox* /*aBox*/, bool /*aStore*/)
    {
      // Computes bounding box.
      std::cout << "ComputeBBox - Not implemented" << std::endl;
    }

    void Extent(UVector3& aMin, UVector3& aMax) const;

  protected: // without description

    // Here are our parameters

    double startPhi;    // Starting phi value (0 < phiStart < 2pi)
    double endPhi;      // end phi value (0 < endPhi-phiStart < 2pi)
    bool   phiIsOpen;  // true if there is a phi segment
    int  numCorner;   // number RZ points
    UPolyconeSideRZ* corners; // corner r,z points

    // Our quick test

    UEnclosingCylinder* enclosingCylinder;

};

#include "UGenericPolycone.icc"

#endif
