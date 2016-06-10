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
// UBox
//
// Class description:
//
//  A simple box defined by half-lengths on the three axis.
//  The center of the box matches the origin of the local reference frame.
//
// 10.06.11 J.Apostolakis, G.Cosmo, A.Gheata
//          Created from original implementation in Geant4 and ROOT
// --------------------------------------------------------------------

#ifndef USOLIDS_UBox
#define USOLIDS_UBox

#ifndef USOLIDS_VUSolid
#include "VUSolid.hh"
#endif

#ifndef USOLIDS_UUtils
#include "UUtils.hh"
#endif

class UBox : public VUSolid
{

  public:
  UBox() : VUSolid(), fDx(0), fDy(0), fDz(0),fCubicVolume(0.), fSurfaceArea(0.) {}
    UBox(const std::string& name, double dx, double dy, double dz);
    virtual ~UBox();

    UBox(const UBox& rhs);
    UBox& operator=(const UBox& rhs);

    // Copy constructor and assignment operator

    void Set(double dx, double dy, double dz);
    void Set(const UVector3& vec);

    // Accessors and modifiers



    inline double GetXHalfLength() const;
    inline double GetYHalfLength() const;
    inline double GetZHalfLength() const;

    void SetXHalfLength(double dx);
    void SetYHalfLength(double dy);
    void SetZHalfLength(double dz);


    // Navigation methods
    EnumInside     Inside(const UVector3& aPoint) const;

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
    bool Normal(const UVector3& aPoint, UVector3& aNormal) const;
//  void Extent ( EAxisType aAxis, double &aMin, double &aMax ) const;
    void Extent(UVector3& aMin, UVector3& aMax) const;
    inline double Capacity();
    inline double SurfaceArea();
    VUSolid* Clone() const;
    UGeometryType GetEntityType() const;
   
    void    ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}

    // Visualisation
    void GetParametersList(int, double* aArray) const
    {
      aArray[0] = GetXHalfLength();
      aArray[1] = GetYHalfLength();
      aArray[2] = GetZHalfLength();
    }

    UVector3 GetPointOnSurface() const;
    UVector3 GetPointOnEdge() const;

    std::ostream& StreamInfo(std::ostream& os) const;


  private:
    double                fDx;   // Half-length on X
    double                fDy;   // Half-length on Y
    double                fDz;   // Half-length on Z
    double       fCubicVolume;   // Cubic Volume
    double       fSurfaceArea;   // Surface Area

};


inline double UBox::GetXHalfLength() const
{
  return fDx;
}
inline double UBox::GetYHalfLength() const
{
  return fDy;
}
inline double UBox::GetZHalfLength() const
{
  return fDz;
}

inline double UBox::Capacity()
{
  if (fCubicVolume != 0.)
  {
    ;
  }
  else
  {
    fCubicVolume = 8 * fDx * fDy * fDz;
  }
  return fCubicVolume;
}

inline double UBox::SurfaceArea()
{
  if (fSurfaceArea != 0.)
  {
    ;
  }
  else
  {
    fSurfaceArea = 8 * (fDx * fDy + fDx * fDz + fDy * fDz);
  }
  return fSurfaceArea;
}


#endif
