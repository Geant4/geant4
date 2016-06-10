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
// 19.09.12 Marek Gayer
//          Created from original implementation in CLHEP
// --------------------------------------------------------------------

#include <cmath>
#include <cstring>

#include "UTransform3D.hh"
#include "UUtils.hh"

const static double kIdRot[9] =
{
  1.0, 0.0, 0.0,
  0.0, 1.0, 0.0,
  0.0, 0.0, 1.0
};

//______________________________________________________________________________
UTransform3D::UTransform3D()
{
  // Dummy constructor
  fTr.Set(0);
  std::memcpy(fRot, kIdRot, sizeof(kIdRot));
}

//______________________________________________________________________________
UTransform3D::UTransform3D(double tx, double ty, double tz,
                           double phi, double theta, double psi)
{
  // Constructor providing a translation and Euler angles
  // See description for SetAngles() method.
  // This represent the composition of : first a rotation about Z axis with
  // angle phi, then a rotation with theta about the rotated X axis, and
  // finally a rotation with psi about the new Z axis.

  fTr.Set(tx, ty, tz);
  SetAngles(phi, theta, psi);
}

//______________________________________________________________________________
UTransform3D::UTransform3D(const UTransform3D& other)
{
  // Copy constructor.
 
    fTr = other.fTr;
    std::memcpy(fRot, other.fRot, sizeof(kIdRot));
 
}

//______________________________________________________________________________
UTransform3D& UTransform3D::operator = (const UTransform3D& other)
{
  if (&other == this) return *this;
  fTr = other.fTr;
  std::memcpy(fRot, other.fRot, sizeof(kIdRot));
  return *this;
}

//______________________________________________________________________________
void UTransform3D::SetAngles(double phi, double theta, double psi)
{
  // Set the rotation from Euler angles in the X-axis convention
  // See: http://mathworld.wolfram.com/EulerAngles.html
  // This represent the composition of : first a rotation about Z axis with
  // angle phi, then a rotation with theta about the rotated X axis, and
  // finally a rotation with psi about the new Z axis.

  // NOTE: angles are in degrees
  double degrad = UUtils::kDegToRad;
  double sinphi = std::sin(degrad * phi);
  double cosphi = std::cos(degrad * phi);
  double sinthe = std::sin(degrad * theta);
  double costhe = std::cos(degrad * theta);
  double sinpsi = std::sin(degrad * psi);
  double cospsi = std::cos(degrad * psi);

  fRot[0] =  cospsi * cosphi - costhe * sinphi * sinpsi;
  fRot[1] = -sinpsi * cosphi - costhe * sinphi * cospsi;
  fRot[2] =  sinthe * sinphi;
  fRot[3] =  cospsi * sinphi + costhe * cosphi * sinpsi;
  fRot[4] = -sinpsi * sinphi + costhe * cosphi * cospsi;
  fRot[5] = -sinthe * cosphi;
  fRot[6] =  sinpsi * sinthe;
  fRot[7] =  cospsi * sinthe;
  fRot[8] =  costhe;
}

//______________________________________________________________________________
void UTransform3D::RotateX(double angle)
{
  // Rotate the transformation about the X axis with a given angle (in degrees).
  double phi = angle * UUtils::kDegToRad;
  double c = std::cos(phi);
  double s = std::sin(phi);
  double v[9];
  v[0] = fRot[0];
  v[1] = fRot[1];
  v[2] = fRot[2];
  v[3] = c * fRot[3] - s * fRot[6];
  v[4] = c * fRot[4] - s * fRot[7];
  v[5] = c * fRot[5] - s * fRot[8];
  v[6] = s * fRot[3] + c * fRot[6];
  v[7] = s * fRot[4] + c * fRot[7];
  v[8] = s * fRot[5] + c * fRot[8];
  std::memcpy(fRot, v, sizeof(kIdRot));

  fTr.Set(fTr.x(), c * fTr.y() - s * fTr.z(), s * fTr.y() + c * fTr.z());
}

//______________________________________________________________________________
void UTransform3D::RotateY(double angle)
{
  // Rotate the transformation about the Y axis with a given angle (in degrees).
  double phi = angle * UUtils::kDegToRad;
  double c = std::cos(phi);
  double s = std::sin(phi);
  double v[9];
  v[0] = c * fRot[0] + s * fRot[6];
  v[1] = c * fRot[1] + s * fRot[7];
  v[2] = c * fRot[2] + s * fRot[8];
  v[3] = fRot[3];
  v[4] = fRot[4];
  v[5] = fRot[5];
  v[6] = -s * fRot[0] + c * fRot[6];
  v[7] = -s * fRot[1] + c * fRot[7];
  v[8] = -s * fRot[2] + c * fRot[8];
  std::memcpy(fRot, v, sizeof(kIdRot));

  fTr.Set(c * fTr.x() + s * fTr.z(), fTr.y(), -s * fTr.x() + c * fTr.z());
}

//______________________________________________________________________________
void UTransform3D::RotateZ(double angle)
{
  // Rotate the transformation about the Z axis with a given angle (in degrees).
  double phi = angle * UUtils::kDegToRad;
  double c = std::cos(phi);
  double s = std::sin(phi);
  double v[9];
  v[0] = c * fRot[0] - s * fRot[3];
  v[1] = c * fRot[1] - s * fRot[4];
  v[2] = c * fRot[2] - s * fRot[5];
  v[3] = s * fRot[0] + c * fRot[3];
  v[4] = s * fRot[1] + c * fRot[4];
  v[5] = s * fRot[2] + c * fRot[5];
  v[6] = fRot[6];
  v[7] = fRot[7];
  v[8] = fRot[8];
  std::memcpy(&fRot[0], v, sizeof(kIdRot));

  fTr.Set(c * fTr.x() - s * fTr.y(), s * fTr.x() + c * fTr.y(), fTr.z());
}

//______________________________________________________________________________
UVector3 UTransform3D::GlobalPoint(const UVector3& local) const
{
  // Returns global point coordinates converted from the local frame defined
  // by the transformation. This is defined by multiplying this transformation
  // with the local vector.
  UVector3 global;
  global.x() = fTr.x() + local.x() * fRot[0] + local.y() * fRot[1] + local.z() * fRot[2];
  global.y() = fTr.y() + local.x() * fRot[3] + local.y() * fRot[4] + local.z() * fRot[5];
  global.z() = fTr.z() + local.x() * fRot[6] + local.y() * fRot[7] + local.z() * fRot[8];
  return global;
}

//______________________________________________________________________________
UVector3 UTransform3D::GlobalVector(const UVector3& local) const
{
  // Returns vector components converted from the local frame defined by the
  // transformation to the global one. This is defined by multiplying this
  // transformation with the local vector while ignoring the translation.
  UVector3 global(
    local.x() * fRot[0] + local.y() * fRot[1] + local.z() * fRot[2],
    local.x() * fRot[3] + local.y() * fRot[4] + local.z() * fRot[5],
    local.x() * fRot[6] + local.y() * fRot[7] + local.z() * fRot[8]);

  return global;
}

//______________________________________________________________________________
UVector3 UTransform3D::LocalPoint(const UVector3& global) const
{
  // Returns local point coordinates converted from the global frame defined
  // by the transformation. This is defined by multiplying the inverse
  // transformation with the global vector.
  UVector3 mt = global - fTr;
  UVector3 local(
    mt.x() * fRot[0] + mt.y() * fRot[3] + mt.z() * fRot[6],
    mt.x() * fRot[1] + mt.y() * fRot[4] + mt.z() * fRot[7],
    mt.x() * fRot[2] + mt.y() * fRot[5] + mt.z() * fRot[8]);
  return local;
}

//______________________________________________________________________________
UVector3 UTransform3D::LocalVector(const UVector3& global) const
{
  // Returns local point coordinates converted from the global frame defined
  // by the transformation. This is defined by multiplying the inverse
  // transformation with the global vector.
  UVector3 local(
    global.x() * fRot[0] + global.y() * fRot[3] + global.z() * fRot[6],
    global.x() * fRot[1] + global.y() * fRot[4] + global.z() * fRot[7],
    global.x() * fRot[2] + global.y() * fRot[5] + global.z() * fRot[8]);

  return local;
}

//______________________________________________________________________________
UTransform3D& UTransform3D::operator *= (const UTransform3D& other)
{
  // Multiply with other transformation.
  fTr.x() = fRot[0] * other.fTr[0] + fRot[1] * other.fTr[1] + fRot[2] * other.fTr[2];
  fTr.y() = fRot[3] * other.fTr[0] + fRot[4] * other.fTr[1] + fRot[5] * other.fTr[2];
  fTr.z() = fRot[6] * other.fTr[0] + fRot[7] * other.fTr[1] + fRot[8] * other.fTr[2];

  double newrot[9];
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      newrot[3 * i + j] = fRot[3 * i] * other.fRot[j] +
                          fRot[3 * i + 1] * other.fRot[3 + j] +
                          fRot[3 * i + 2] * other.fRot[6 + j];
    }
  }
  std::memcpy(fRot, newrot, sizeof(kIdRot));
  return *this;
}

//______________________________________________________________________________
UTransform3D& UTransform3D::operator *= (const UVector3& vect)
{
  // Multiply with a vector.
  fTr.x() = fRot[0] * vect.x() + fRot[1] * vect.y() + fRot[2] * vect.z();
  fTr.y() = fRot[3] * vect.x() + fRot[4] * vect.y() + fRot[5] * vect.z();
  fTr.z() = fRot[6] * vect.x() + fRot[7] * vect.y() + fRot[8] * vect.z();
  return *this;
}

//______________________________________________________________________________
UVector3 operator * (const UVector3& /*p*/, const UTransform3D& /*trans*/)
{
  // Multiply matrix with translation.
  UVector3 vect;
  return vect;
}
