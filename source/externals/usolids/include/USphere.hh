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
// USphere
//
// Class description:
//
//   A USphere is, in the general case, a section of a spherical shell,
//   between specified phi and theta angles
//
//   The phi and theta segments are described by a starting angle,
//   and the +ve delta angle for the shape.
//   If the delta angle is >=2*UUtils::kPi, or >=UUtils::kPi the shape is treated as
//   continuous in phi or theta respectively.
//
//   Theta must lie between 0-UUtils::kPi (incl).
//
//   Member Data:
//
//   fRmin  inner radius
//   fRmax  outer radius
//
//   fSPhi  starting angle of the segment in radians
//   fDPhi  delta angle of the segment in radians
//
//   fSTheta  starting angle of the segment in radians
//   fDTheta  delta angle of the segment in radians
//
//
//   Note:
//      Internally fSPhi & fDPhi are adjusted so that fDPhi<=2PI,
//      and fDPhi+fSPhi<=2PI. This enables simpler comparisons to be
//      made with (say) Phi of a point.
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef USphere_HH
#define USphere_HH

#include <sstream>

#include "VUSolid.hh"

class UVisExtent;

class USphere : public VUSolid
{
  public: // with description

    USphere(const std::string& pName,
            double pRmin, double pRmax,
            double pSPhi, double pDPhi,
            double pSTheta, double pDTheta);
    //
    // Constructs a sphere or sphere shell section
    // with the given name and dimensions

    ~USphere();
    //
    // Destructor

    // Accessors

    inline double GetInnerRadius() const;
    inline double GetOuterRadius() const;
    inline double GetStartPhiAngle() const;
    inline double GetDeltaPhiAngle() const;
    inline double GetStartThetaAngle() const;
    inline double GetDeltaThetaAngle() const;

    // Modifiers

    inline void SetInnerRadius(double newRMin);
    inline void SetOuterRadius(double newRmax);
    inline void SetStartPhiAngle(double newSphi, bool trig = true);
    inline void SetDeltaPhiAngle(double newDphi);
    inline void SetStartThetaAngle(double newSTheta);
    inline void SetDeltaThetaAngle(double newDTheta);

    // Methods for solid

    inline double Capacity();
    double SurfaceArea();


    VUSolid::EnumInside Inside(const UVector3& p) const;

    bool Normal(const UVector3& p, UVector3& n) const;

    double DistanceToIn(const UVector3& p, const UVector3& v, double aPstep = UUtils::kInfinity) const;

    double SafetyFromOutside(const UVector3& p, bool aAccurate = false) const;

    double DistanceToOut(const UVector3& p, const UVector3& v, UVector3& n, bool& validNorm, double aPstep = UUtils::kInfinity) const;


    double SafetyFromInside(const UVector3& p, bool aAccurate = false) const;

    UGeometryType GetEntityType() const;

    UVector3 GetPointOnSurface() const;

    VUSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

    // Visualisation functions

    UVisExtent   GetExtent() const;



  public: // without description

    void Extent(UVector3& aMin, UVector3& aMax) const;

    void GetParametersList(int /*aNumber*/, double* /*aArray*/) const;


    virtual void    ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}

    USphere(const USphere& rhs);
    USphere& operator=(const USphere& rhs);
    // Copy constructor and assignment operator.

    // Old access functions

    inline double GetRmin()  const;
    inline double GetRmax()  const;
    inline double GetSPhi()  const;
    inline double GetDPhi()  const;
    inline double GetSTheta() const;
    inline double GetDTheta() const;
    inline double GetInsideRadius() const;
    inline void SetInsideRadius(double newRmin);

  private:

    double fCubicVolume;
    double fSurfaceArea;
    inline void Initialize();
    //
    // Reset relevant values to zero

    inline void CheckThetaAngles(double sTheta, double dTheta);
    inline void CheckSPhiAngle(double sPhi);
    inline void CheckDPhiAngle(double dPhi);
    inline void CheckPhiAngles(double sPhi, double dPhi);
    //
    // Reset relevant flags and angle values

    inline void InitializePhiTrigonometry();
    inline void InitializeThetaTrigonometry();
    //
    // Recompute relevant trigonometric values and cache them

    UVector3 ApproxSurfaceNormal(const UVector3& p) const;
    //
    // Algorithm for SurfaceNormal() following the original
    // specification for points not on the surface

  private:

    // Used by distanceToOut
    //
    enum ESide {kNull, kRMin, kRMax, kSPhi, kEPhi, kSTheta, kETheta};

    // used by normal
    //
    enum ENorm {kNRMin, kNRMax, kNSPhi, kNEPhi, kNSTheta, kNETheta};

    double fRminTolerance, kTolerance, kAngTolerance,
           kRadTolerance, fEpsilon;
    //
    // Radial and angular tolerances

    double fRmin, fRmax, fSPhi, fDPhi, fSTheta, fDTheta;
    //
    // Radial and angular dimensions

    double sinCPhi, cosCPhi, cosHDPhiOT, cosHDPhiIT,
           sinSPhi, cosSPhi, sinEPhi, cosEPhi, hDPhi, cPhi, ePhi;
    //
    // Cached trigonometric values for Phi angle

    double sinSTheta, cosSTheta, sinETheta, cosETheta,
           tanSTheta, tanSTheta2, tanETheta, tanETheta2, eTheta;
    //
    // Cached trigonometric values for Theta angle

    bool fFullPhiSphere, fFullThetaSphere, fFullSphere;
    //
    // Flags for identification of section, shell or full sphere
};

inline
double USphere::GetInsideRadius() const
{
  return fRmin;
}

inline
double USphere::GetInnerRadius() const
{
  return fRmin;
}

inline
double USphere::GetOuterRadius() const
{
  return fRmax;
}

inline
double USphere::GetStartPhiAngle() const
{
  return fSPhi;
}

inline
double USphere::GetDeltaPhiAngle() const
{
  return fDPhi;
}

inline
double USphere::GetStartThetaAngle() const
{
  return fSTheta;
}

double USphere::GetDeltaThetaAngle() const
{
  return fDTheta;
}

inline
void USphere::Initialize()
{
  fCubicVolume = 0.;
  fSurfaceArea = 0.;
}

inline
void USphere::InitializePhiTrigonometry()
{
  hDPhi = 0.5 * fDPhi;        // half delta phi
  cPhi  = fSPhi + hDPhi;
  ePhi  = fSPhi + fDPhi;

  sinCPhi   = std::sin(cPhi);
  cosCPhi   = std::cos(cPhi);
  cosHDPhiIT = std::cos(hDPhi - 0.5 * kAngTolerance); // inner/outer tol half dphi
  cosHDPhiOT = std::cos(hDPhi + 0.5 * kAngTolerance);
  sinSPhi = std::sin(fSPhi);
  cosSPhi = std::cos(fSPhi);
  sinEPhi = std::sin(ePhi);
  cosEPhi = std::cos(ePhi);
}

inline
void USphere::InitializeThetaTrigonometry()
{
  eTheta  = fSTheta + fDTheta;

  sinSTheta = std::sin(fSTheta);
  cosSTheta = std::cos(fSTheta);
  sinETheta = std::sin(eTheta);
  cosETheta = std::cos(eTheta);

  tanSTheta = std::tan(fSTheta);
  tanSTheta2 = tanSTheta * tanSTheta;
  tanETheta = std::tan(eTheta);
  tanETheta2 = tanETheta * tanETheta;
}

inline
void USphere::CheckThetaAngles(double sTheta, double dTheta)
{
  if ((sTheta < 0) || (sTheta > UUtils::kPi))
  {
    std::ostringstream message;
    message << "sTheta outside 0-PI range." << std::endl
            << "Invalid starting Theta angle for solid: " << GetName();
    UUtils::Exception("USphere::CheckThetaAngles()", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }
  else
  {
    fSTheta = sTheta;
  }
  if (dTheta + sTheta >= UUtils::kPi)
  {
    fDTheta = UUtils::kPi - sTheta;
  }
  else if (dTheta > 0)
  {
    fDTheta = dTheta;
  }
  else
  {
    std::ostringstream message;
    message << "Invalid dTheta." << std::endl
            << "Negative delta-Theta (" << dTheta << "), for solid: "
            << GetName();
    UUtils::Exception("USphere::CheckThetaAngles()", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }
  if (fDTheta - fSTheta < UUtils::kPi)
  {
    fFullThetaSphere = false;
  }
  else
  {
    fFullThetaSphere = true ;
  }
  fFullSphere = fFullPhiSphere && fFullThetaSphere;

  InitializeThetaTrigonometry();
}

inline
void USphere::CheckSPhiAngle(double sPhi)
{
  // Ensure fSphi in 0-2PI or -2PI-0 range if shape crosses 0

  if (sPhi < 0)
  {
    fSPhi = 2 * UUtils::kPi - std::fmod(std::fabs(sPhi), 2 * UUtils::kPi);
  }
  else
  {
    fSPhi = std::fmod(sPhi, 2 * UUtils::kPi) ;
  }
  if (fSPhi + fDPhi > 2 * UUtils::kPi)
  {
    fSPhi -= 2 * UUtils::kPi ;
  }
}

inline
void USphere::CheckDPhiAngle(double dPhi)
{
  fFullPhiSphere = true;
  if (dPhi >= 2 * UUtils::kPi - kAngTolerance * 0.5)
  {
    fDPhi = 2 * UUtils::kPi;
    fSPhi = 0;
  }
  else
  {
    fFullPhiSphere = false;
    if (dPhi > 0)
    {
      fDPhi = dPhi;
    }
    else
    {
      std::ostringstream message;
      message << "Invalid dphi." << std::endl
              << "Negative delta-Phi (" << dPhi << "), for solid: "
              << GetName();
      UUtils::Exception("USphere::CheckDPhiAngle()", "GeomSolids0002",
                        UFatalError, 1, message.str().c_str());
    }
  }
}

inline
void USphere::CheckPhiAngles(double sPhi, double dPhi)
{
  CheckDPhiAngle(dPhi);
  //if (!fFullPhiSphere && sPhi) { CheckSPhiAngle(sPhi); }
  if (!fFullPhiSphere)
  {
    CheckSPhiAngle(sPhi);
  }
  fFullSphere = fFullPhiSphere && fFullThetaSphere;

  InitializePhiTrigonometry();
}

inline
void USphere::SetInsideRadius(double newRmin)
{
  fRmin = newRmin;
  fRminTolerance = (fRmin) ? std::max(kRadTolerance, fEpsilon * fRmin) : 0;
  Initialize();
}

inline
void USphere::SetInnerRadius(double newRmin)
{
  SetInsideRadius(newRmin);
}

inline
void USphere::SetOuterRadius(double newRmax)
{
  fRmax = newRmax;
  kTolerance = std::max(kRadTolerance, fEpsilon * fRmax);
  Initialize();
}

inline
void USphere::SetStartPhiAngle(double newSPhi, bool compute)
{
  // Flag 'compute' can be used to explicitely avoid recomputation of
  // trigonometry in case SetDeltaPhiAngle() is invoked afterwards

  CheckSPhiAngle(newSPhi);
  fFullPhiSphere = false;
  if (compute)
  {
    InitializePhiTrigonometry();
  }
  Initialize();
}

inline
void USphere::SetDeltaPhiAngle(double newDPhi)
{
  CheckPhiAngles(fSPhi, newDPhi);
  Initialize();
}

inline
void USphere::SetStartThetaAngle(double newSTheta)
{
  CheckThetaAngles(newSTheta, fDTheta);
  Initialize();
}

inline
void USphere::SetDeltaThetaAngle(double newDTheta)
{
  CheckThetaAngles(fSTheta, newDTheta);
  Initialize();
}

// Old access functions

inline
double USphere::GetRmin() const
{
  return GetInsideRadius();
}

inline
double USphere::GetRmax() const
{
  return GetOuterRadius();
}

inline
double USphere::GetSPhi() const
{
  return GetStartPhiAngle();
}

inline
double USphere::GetDPhi() const
{
  return GetDeltaPhiAngle();
}

inline
double USphere::GetSTheta() const
{
  return GetStartThetaAngle();
}

inline
double USphere::GetDTheta() const
{
  return GetDeltaThetaAngle();
}

inline
double USphere::Capacity()
{
  if (fCubicVolume != 0.)
  {
    ;
  }
  else
  {
    fCubicVolume = fDPhi * (std::cos(fSTheta) - std::cos(fSTheta + fDTheta)) *
                   (fRmax * fRmax * fRmax - fRmin * fRmin * fRmin) / 3.;
  }
  return fCubicVolume;
}

#endif
