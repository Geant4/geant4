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
// UTrd
//
// Class description:
//
//   A UTrd is a trapezoid with the x and y dimensions varying along z
//   functions.
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef USOLIDS_UTrd
#define USOLIDS_UTrd

#include "VUSolid.hh"
#include "UUtils.hh"

class UTrd : public VUSolid
{
    enum ESide {kUndefined, kPX, kMX, kPY, kMY, kPZ, kMZ};
  public:
    UTrd() : VUSolid(), fDx1(0), fDx2(0), fDy1(0), fDy2(0), fDz(0) {}
    UTrd(const std::string& pName, double pdx1, double pdx2, double pdy1, double pdy2, double pdz);
    virtual ~UTrd() {}

    UTrd(const UTrd& rhs);
    UTrd& operator=(const UTrd& rhs);

    // Copy constructor and assignment operator

    // Accessors

    inline double GetXHalfLength1() const;
    inline double GetXHalfLength2() const;
    inline double GetYHalfLength1() const;
    inline double GetYHalfLength2() const;
    inline double GetZHalfLength()  const;

    // Modifiers

    inline void SetXHalfLength1(double val);
    inline void SetXHalfLength2(double val);
    inline void SetYHalfLength1(double val);
    inline void SetYHalfLength2(double val);
    inline void SetZHalfLength(double val);
    // Navigation methods
    EnumInside     Inside(const UVector3& aPoint) const;

    virtual double SafetyFromInside(const UVector3& aPoint, bool aAccurate = false) const;

    double SafetyFromInsideAccurate(const UVector3& aPoint) const;

    virtual double SafetyFromOutside(const UVector3& aPoint, bool aAccurate = false) const;

    double SafetyFromOutsideAccurate(const UVector3& aPoint) const;

    virtual double DistanceToIn(const UVector3& aPoint,
                                const UVector3& aDirection,
                                // UVector3       &aNormalVector,
                                double aPstep = UUtils::kInfinity) const;


    virtual double DistanceToOut(const UVector3& aPoint,
                                 const UVector3& aDirection,
                                 UVector3&       aNormalVector,
                                 bool&           aConvex,
                                 double aPstep = UUtils::kInfinity) const;

    virtual bool Normal(const UVector3& aPoint, UVector3& aNormal) const;

    void CheckAndSetAllParameters ( double pdx1, double pdx2,
                                    double pdy1, double pdy2,
                                    double pdz );

    void SetAllParameters ( double pdx1, double pdx2,
                            double pdy1, double pdy2,
                            double pdz );

//  virtual void Extent ( EAxisType aAxis, double &aMin, double &aMax ) const;
    void Extent(UVector3& aMin, UVector3& aMax) const;
    inline double Capacity();
    inline double SurfaceArea();
    VUSolid* Clone() const;
    UGeometryType GetEntityType() const;
 
    virtual void ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}

    //G4Visualisation
    virtual void GetParametersList(int /*aNumber*/, double* /*aArray*/) const;
    std::ostream& StreamInfo(std::ostream& os) const;

    UVector3 GetPointOnSurface() const;

  private:
    UVector3 ApproxSurfaceNormal(const UVector3& p) const;
    inline double amin(int n, const double* a) const;
    inline double amax(int n, const double* a)const;
    double fDx1, fDx2, fDy1, fDy2, fDz;
    double       fCubicVolume;   // Cubic Volume
    double       fSurfaceArea;   // Surface Area
};

#include "UTrd.icc"

#endif
