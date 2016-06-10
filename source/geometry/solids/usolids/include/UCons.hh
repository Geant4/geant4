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
// UCons
//
// Class description:
//
//   A UCons is, in the general case, a Phi segment of a cone, with
//   half-length fDz, inner and outer radii specified at -fDz and +fDz.
//   The Phi segment is described by a starting fSPhi angle, and the
//   +fDPhi delta angle for the shape.
//   If the delta angle is >=2*UUtils::kPi, the shape is treated as
//   continuous in Phi
//
//   Member Data:
//
//   fRmin1  inside radius at  -fDz
//   fRmin2  inside radius at  +fDz
//   fRmax1  outside radius at -fDz
//   fRmax2  outside radius at +fDz
//   fDz half length in z
//
//   fSPhi starting angle of the segment in radians
//   fDPhi delta angle of the segment in radians
//
//   fPhiFullCone   Boolean variable used for indicate the Phi Section
//
//   Note:
//      Internally fSPhi & fDPhi are adjusted so that fDPhi<=2PI,
//      and fDPhi+fSPhi<=2PI. This enables simpler comparisons to be
//      made with (say) Phi of a point.
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UCons_HH
#define UCons_HH

#include "VUSolid.hh"

class UCons : public VUSolid
{
  public: // with description

    UCons(const std::string& pName,
          double pRmin1, double pRmax1,
          double pRmin2, double pRmax2,
          double pDz,
          double pSPhi, double pDPhi);
    //
    // Constructs a cone with the given name and dimensions

    ~UCons() ;
    //
    // Destructor

    // Accessors

    inline double GetInnerRadiusMinusZ() const;
    inline double GetOuterRadiusMinusZ() const;
    inline double GetInnerRadiusPlusZ() const;
    inline double GetOuterRadiusPlusZ() const;
    inline double GetZHalfLength()       const;
    inline double GetStartPhiAngle()     const;
    inline double GetDeltaPhiAngle()     const;

    // Modifiers

    inline void SetInnerRadiusMinusZ(double Rmin1);
    inline void SetOuterRadiusMinusZ(double Rmax1);
    inline void SetInnerRadiusPlusZ(double Rmin2);
    inline void SetOuterRadiusPlusZ(double Rmax2);
    inline void SetZHalfLength(double newDz);
    inline void SetStartPhiAngle(double newSPhi, bool trig = true);
    inline void SetDeltaPhiAngle(double newDPhi);

    // Other methods for solid

    inline double GetCubicVolume();
    inline double GetSurfaceArea();


//    inline VUSolid::EnumInside Inside( const UVector3& p ) const;

    bool Normal(const UVector3& p, UVector3& n) const;

    double DistanceToIn(const UVector3& p, const UVector3& v, double aPstep = UUtils::kInfinity) const;

    double SafetyFromOutside(const UVector3& p, bool precise) const;



    double DistanceToOut(const UVector3& aPoint,
                         const UVector3&  aDirection,
                         UVector3&       aNormalVector,
                         bool&           aConvex,
                         double aPstep = UUtils::kInfinity) const;

    double SafetyFromInside(const UVector3& p, bool precise) const;

    UGeometryType GetEntityType() const;

    UVector3 GetPointOnSurface() const;

    VUSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

//      void Extent (EAxisType aAxis, double &aMin, double &aMax) const;
    void Extent(UVector3& aMin, UVector3& aMax) const;

    virtual void GetParametersList(int /*aNumber*/, double* /*aArray*/) const;

    virtual void ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}

    // Visualisation functions


    inline VUSolid::EnumInside Inside(const UVector3& p) const
    {
      double r2, rl, rh, pPhi, tolRMin, tolRMax; // rh2, rl2;
      VUSolid::EnumInside in;
      static const double halfCarTolerance = VUSolid::Tolerance() * 0.5;
      static const double halfRadTolerance = kRadTolerance * 0.5;
      static const double halfAngTolerance = kAngTolerance * 0.5;

      if (std::fabs(p.z) > fDz + halfCarTolerance)
      {
        return in = eOutside;
      }
      else if (std::fabs(p.z) >= fDz - halfCarTolerance)
      {
        in = eSurface;
      }
      else
      {
        in = eInside;
      }

      r2 = p.x * p.x + p.y * p.y;
      rl = 0.5 * (fRmin2 * (p.z + fDz) + fRmin1 * (fDz - p.z)) / fDz;
      rh = 0.5 * (fRmax2 * (p.z + fDz) + fRmax1 * (fDz - p.z)) / fDz;

      // rh2 = rh*rh;

      tolRMin = rl - halfRadTolerance;
      if (tolRMin < 0)
      {
        tolRMin = 0;
      }
      tolRMax = rh + halfRadTolerance;

      if ((r2 < tolRMin * tolRMin) || (r2 > tolRMax * tolRMax))
      {
        return in = eOutside;
      }

      if (rl)
      {
        tolRMin = rl + halfRadTolerance;
      }
      else
      {
        tolRMin = 0.0;
      }
      tolRMax = rh - halfRadTolerance;

      if (in == eInside) // else it's eSurface already
      {
        if ((r2 < tolRMin * tolRMin) || (r2 >= tolRMax * tolRMax))
        {
          in = eSurface;
        }
      }
      if (!fPhiFullCone && ((p.x != 0.0) || (p.y != 0.0)))
      {
        pPhi = std::atan2(p.y, p.x);

        if (pPhi < fSPhi - halfAngTolerance)
        {
          pPhi += 2 * UUtils::kPi;
        }
        else if (pPhi > fSPhi + fDPhi + halfAngTolerance)
        {
          pPhi -= 2 * UUtils::kPi;
        }

        if ((pPhi < fSPhi - halfAngTolerance) ||
            (pPhi > fSPhi + fDPhi + halfAngTolerance))
        {
          return in = eOutside;
        }

        else if (in == eInside) // else it's eSurface anyway already
        {
          if ((pPhi < fSPhi + halfAngTolerance) ||
              (pPhi > fSPhi + fDPhi - halfAngTolerance))
          {
            in = eSurface;
          }
        }
      }
      else if (!fPhiFullCone)
      {
        in = eSurface;
      }

      return in;
    }

  public: // without description

    UCons();
    //
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

    UCons(const UCons& rhs);
    UCons& operator=(const UCons& rhs);
    // Copy constructor and assignment operator.

    //  Old access functions

    inline double   GetRmin1() const;
    inline double   GetRmax1() const;
    inline double   GetRmin2() const;
    inline double   GetRmax2() const;
    inline double   GetDz()   const;
    inline double   GetSPhi() const;
    inline double   GetDPhi() const;

  private:

    double fCubicVolume, fSurfaceArea;

    inline double Capacity();
    inline double SurfaceArea();

    inline void Initialize();
    //
    // Reset relevant values to zero

    inline void CheckSPhiAngle(double sPhi);
    inline void CheckDPhiAngle(double dPhi);
    inline void CheckPhiAngles(double sPhi, double dPhi);
    //
    // Reset relevant flags and angle values

    inline void InitializeTrigonometry();
    //
    // Recompute relevant trigonometric values and cache them

    UVector3 ApproxSurfaceNormal(const UVector3& p) const;
    //
    // Algorithm for SurfaceNormal() following the original
    // specification for points not on the surface

  private:

    // Used by distanceToOut
    //
    enum ESide {kNull, kRMin, kRMax, kSPhi, kEPhi, kPZ, kMZ};

    // used by normal
    //
    enum ENorm {kNRMin, kNRMax, kNSPhi, kNEPhi, kNZ};

    double kRadTolerance, kAngTolerance;
    //
    // Radial and angular tolerances

    double fRmin1, fRmin2, fRmax1, fRmax2, fDz, fSPhi, fDPhi;
    //
    // Radial and angular dimensions

    double sinCPhi, cosCPhi, cosHDPhiOT, cosHDPhiIT,
           sinSPhi, cosSPhi, sinEPhi, cosEPhi;
    //
    // Cached trigonometric values

    bool fPhiFullCone;

    double secRMin, tanRMin, tanRMax, secRMax;

//    double fSinPhi;

    //
    // Flag for identification of section or full cone
};

#include "UCons.icc"

#endif
