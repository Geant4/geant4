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
// UTet
//
// Class description:
//
//   A UTet is a tetrahedrasolid.
//
// 19.07.13 Tatiana Nikitina
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UTet_hh
#define UTet_hh

#include "VUSolid.hh"
#include "UUtils.hh"

class UTet : public VUSolid
{

  public:  // with description

    UTet(const std::string& name,
         UVector3 anchor,
         UVector3 p2,
         UVector3 p3,
         UVector3 p4,
         bool* degeneracyFlag = 0);

    virtual ~UTet();


    // Methods for solid

    EnumInside Inside(const UVector3&  p) const;

    bool Normal(const UVector3& aPoint, UVector3& aNormal) const;

    double  SafetyFromInside(const UVector3& aPoint,
                             bool aAccurate = false) const;
    double  SafetyFromOutside(const UVector3& aPoint,
                              bool aAccurate = false) const;
    double  DistanceToIn(const UVector3& aPoint,
                         const UVector3& aDirection,
                         // UVector3       &aNormalVector,
                         double aPstep = UUtils::kInfinity) const;

    double DistanceToOut(const UVector3& aPoint,
                         const UVector3& aDirection,
                         UVector3&       aNormalVector,
                         bool&           aConvex,
                         double aPstep = UUtils::kInfinity) const;

    void Extent(UVector3& aMin, UVector3& aMax) const;
    
    double Capacity();
    double SurfaceArea();
    
    UGeometryType GetEntityType() const;
  
    void ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}

    // Visualisation
    void GetParametersList(int aNumber, double* aArray) const;

    VUSolid* Clone() const;

    UVector3 GetPointOnSurface() const;
    std::ostream& StreamInfo(std::ostream& os) const;

  public:   // without description

    UTet(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

    UTet(const UTet& rhs);
    UTet& operator=(const UTet& rhs);
    // Copy constructor and assignment operator.

    void PrintWarnings(bool flag)
    {
      warningFlag = flag;
    }
    static bool CheckDegeneracy(UVector3& anchor,
                                UVector3& p2,
                                UVector3& p3,
                                UVector3& p4);
    std::vector<UVector3> GetVertices() const;
    // Return the four vertices of the shape.

  private:

    double fCubicVolume, fSurfaceArea;

    UVector3 GetPointOnFace(UVector3 p1, UVector3 p2,
                            UVector3 p3, double& area) const;
    static const char CVSVers[];

  private:

    UVector3 fAnchor, fP2, fP3, fP4, fMiddle;
    UVector3 fNormal123, fNormal142, fNormal134, fNormal234;

    bool warningFlag;

    double fCdotN123, fCdotN142, fCdotN134, fCdotN234;
    double fXMin, fXMax, fYMin, fYMax, fZMin, fZMax;
    double fDx, fDy, fDz, fTol, fMaxSize;
};

#endif
