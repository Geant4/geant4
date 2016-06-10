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
// UTransform3D
//
// Class description:
//
//  UTransform3D: General transformation made by rotation + translation
//
// 19.10.12 Marek Gayer
//          Created from original implementation in CLHEP
// --------------------------------------------------------------------

#ifndef USOLIDS_UTransform3D
#define USOLIDS_UTransform3D

#include "UVector3.hh"

class UTransform3D
{
  public:

    UTransform3D();  // Initialize to identity
    UTransform3D(double tx, double ty, double tz,
                 double phi = 0., double theta = 0., double psi = 0.);
    UTransform3D(const UTransform3D& other);
    ~UTransform3D() {}

    void                 RotateX(double angle);
    void                 RotateY(double angle);
    void                 RotateZ(double angle);
    void                 SetAngles(double phi, double theta, double psi);

    // Local<->global coordinate and vector conversions
    UVector3             GlobalPoint(const UVector3& local) const;
    UVector3             GlobalVector(const UVector3& local) const;
    UVector3             LocalPoint(const UVector3& global) const;
    UVector3             LocalVector(const UVector3& global) const;

    // Operators
    UTransform3D& operator = (const UTransform3D& other);
    UTransform3D& operator *= (const UTransform3D& other);
    UTransform3D& operator *= (const UVector3& vect);

    UVector3          fTr;       // Translation
    double            fRot[9];   // Rotation
};
// Vector-matrix multiplication
UVector3 operator * (const UVector3& p, const UTransform3D& trans);

#endif
