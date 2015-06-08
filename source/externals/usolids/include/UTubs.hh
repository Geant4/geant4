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
// UTubs
//
// Class description:
//
//   A tube or tube segment with curved sides parallel to
//   the z-axis. The tube has a specified half-length along
//   the z-axis, about which it is centered, and a given
//   minimum and maximum radius. A minimum radius of 0
//   corresponds to filled tube /cylinder. The tube segment is
//   specified by starting and delta angles for phi, with 0
//   being the +x axis, PI/2 the +y axis.
//   A delta angle of 2PI signifies a complete, unsegmented
//   tube/cylinder.
//
//   Member Data:
//
//   fRMin  Inner radius
//   fRMax  Outer radius
//   fDz  half length in z
//
//   fSPhi  The starting phi angle in radians,
//          adjusted such that fSPhi+fDPhi<=2PI, fSPhi>-2PI
//
//   fDPhi  Delta angle of the segment.
//
//   fPhiFullTube  Boolean variable used for indicate the Phi Section
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UTUBS_HH
#define UTUBS_HH

#include "VUSolid.hh"
#include <sstream>

class UTubs : public VUSolid
{
  public:

    UTubs(const std::string& pName,
          double pRMin,
          double pRMax,
          double pDz,
          double pSPhi,
          double pDPhi);
    //
    // Constructs a tubs with the given name and dimensions

    virtual ~UTubs();
    //
    // Destructor

    // Accessors

    inline double GetInnerRadius() const;
    inline double GetOuterRadius() const;
    inline double GetZHalfLength() const;
    inline double GetStartPhiAngle() const;
    inline double GetDeltaPhiAngle() const;

    // Modifiers

    inline void SetInnerRadius(double newRMin);
    inline void SetOuterRadius(double newRMax);
    inline void SetZHalfLength(double newDz);
    inline void SetStartPhiAngle(double newSPhi, bool trig = true);
    inline void SetDeltaPhiAngle(double newDPhi);

    // Methods for solid

    inline double Capacity();
    inline double SurfaceArea();

    VUSolid::EnumInside Inside(const UVector3& p) const;

    bool Normal(const UVector3& p, UVector3& normal) const;

    double DistanceToIn(const UVector3& p, const UVector3& v,
                        double aPstep = UUtils::kInfinity) const;
    double SafetyFromInside(const UVector3& p, bool precise = false) const;
    double DistanceToOut(const UVector3& p, const UVector3& v, UVector3& n,
                        bool& validNorm, double aPstep=UUtils::kInfinity) const;
    double SafetyFromOutside(const UVector3& p, bool precise = false ) const;
 
    inline double SafetyFromInsideR(const UVector3& p, const double rho,
                                    bool precise = false) const;
    inline double SafetyFromOutsideR(const UVector3& p, const double rho,
                                     bool precise = false) const;
    UGeometryType GetEntityType() const;

    UVector3 GetPointOnSurface() const;

    VUSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

    void Extent(UVector3& aMin, UVector3& aMax) const;

    virtual void GetParametersList(int /*aNumber*/, double* /*aArray*/) const;
    virtual void ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}

  public:

    UTubs();
    //
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

    UTubs(const UTubs& rhs);
    UTubs& operator=(const UTubs& rhs);
    // Copy constructor and assignment operator.

    inline double GetRMin() const;
    inline double GetRMax() const;
    inline double GetDz() const;
    inline double GetSPhi() const;
    inline double GetDPhi() const;

  protected:

    //    UVector3List*
    //    CreateRotatedVertices( const UAffineTransform& pTransform ) const;
    //
    // Creates the List of transformed vertices in the format required
    // for VUSolid:: ClipCrossSection and ClipBetweenSections

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

    virtual UVector3 ApproxSurfaceNormal(const UVector3& p) const;
    //
    // Algorithm for SurfaceNormal() following the original
    // specification for points not on the surface

     inline double SafetyToPhi(const UVector3& p, const double rho, bool& outside) const;

  protected:

    double fCubicVolume, fSurfaceArea;
    // Used by distanceToOut
    //
    enum ESide {kNull, kRMin, kRMax, kSPhi, kEPhi, kPZ, kMZ};

    // Used by normal
    //
    enum ENorm {kNRMin, kNRMax, kNSPhi, kNEPhi, kNZ};

    double kRadTolerance, kAngTolerance;
    //
    // Radial and angular tolerances

    double fRMin, fRMax, fDz, fSPhi, fDPhi;
    //
    // Radial and angular dimensions

    double fSinCPhi, fCosCPhi, fCosHDPhiOT, fCosHDPhiIT,
           fSinSPhi, fCosSPhi, fSinEPhi, fCosEPhi, fSinSPhiDPhi, fCosSPhiDPhi;
    //
    // Cached trigonometric values

    bool fPhiFullTube;
    //
    // Flag for identification of section or full tube
};

#include "UTubs.icc"

#endif
