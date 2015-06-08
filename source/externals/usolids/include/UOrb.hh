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
// UOrb
//
// Class description:
//
//   A simple Orb defined by half-lengths on the three axis.
//   The center of the Orb matches the origin of the local reference frame.
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef USOLIDS_UOrb
#define USOLIDS_UOrb

#include "VUSolid.hh"
#include "UUtils.hh"

class UOrb : public VUSolid
{

  public:
    UOrb() : VUSolid(), fR(0), fRTolerance(0) {}
    UOrb(const std::string& name, double pRmax);
    ~UOrb() {}

    UOrb(const UOrb& rhs);
    UOrb& operator=(const UOrb& rhs);

    // Accessors
    inline double GetRadius() const;
    // Modifiers
    inline void SetRadius(double newRmax);

    // Navigation methods
    EnumInside     Inside(const UVector3& aPo6int) const;

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

    bool Normal(const UVector3& aPoint, UVector3& aNormal) const;
    void Extent(UVector3& aMin, UVector3& aMax) const;
    inline double Capacity();
    inline double SurfaceArea();
    UGeometryType GetEntityType() const;
   
    void ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}

    // Visualisation
    void GetParametersList(int /*aNumber*/, double* /*aArray*/) const;

    VUSolid* Clone() const;

    double GetRadialTolerance()
    {
      return fRTolerance;
    }

    UVector3 GetPointOnSurface() const;

    std::ostream& StreamInfo(std::ostream& os) const;

  private:
    double fR;
    double fRTolerance;
    double fCubicVolume;   // Cubic Volume
    double fSurfaceArea;   // Surface Area

    double DistanceToOutForOutsidePoints(const UVector3& p, const UVector3& v, UVector3& n) const;

};

inline double UOrb::GetRadius() const
{
  return fR;
}
inline void UOrb::SetRadius(double newRmax)
{
  fR = newRmax;
  fCubicVolume = 0.;
  fSurfaceArea = 0.;
}

inline double UOrb::Capacity()
{
  if (fCubicVolume != 0.)
  {
    ;
  }
  else
  {
    fCubicVolume = (4 * UUtils::kPi / 3) * fR * fR * fR;
  }
  return fCubicVolume;
}

inline double UOrb::SurfaceArea()
{
  if (fSurfaceArea != 0.)
  {
    ;
  }
  else
  {
    fSurfaceArea = (4 * UUtils::kPi) * fR * fR;
  }
  return fSurfaceArea;
}

#endif
