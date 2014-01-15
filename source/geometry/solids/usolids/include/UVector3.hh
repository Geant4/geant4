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
// UVector3
//
// Class description:
//
//    Bucket type for Vector type.
//
// 19.09.12 Marek Gayer
//          Created from original implementation in CLHEP
// --------------------------------------------------------------------

#ifndef USOLIDS_UVector3
#define USOLIDS_UVector3

#include <cmath>
#include <iostream>
#include <fstream>

struct UVector3
{
  public:
    UVector3()
    {
      x = y = z = 0.0;
    }
    UVector3(double xval, double yval, double zval)
    {
      x = xval;
      y = yval;
      z = zval;
    }
    UVector3(double theta, double phi);
    UVector3(const double coord[3])
    {
      x = coord[0];
      y = coord[1];
      z = coord[2];
    }

    inline UVector3& operator = (const UVector3& v);
    inline UVector3& operator = (const double* vect);
    // Assignments

    inline bool operator == (const UVector3&) const;
    inline bool operator != (const UVector3&) const;
    // Comparisons.

    inline UVector3 operator - () const;
    // Unary minus.

    inline UVector3& operator += (const UVector3&);
    // Addition.

    inline UVector3& operator -= (const UVector3&);
    // Subtraction.

    inline double& operator[](int index);

    inline double operator[](int index) const;

    inline UVector3& operator *= (double);
    // Scaling with real numbers.

    inline UVector3& operator /= (double);
    // Dividing with real numbers.

    inline double Dot(const UVector3&) const;
    // Scalar product.

    inline UVector3 Cross(const UVector3&) const;
    // Cross product.

    double Angle(const UVector3&) const;
    // The angle w.r.t. another 3-vector.

    UVector3 Unit() const;
    // Unit vector parallel to this.

    inline bool IsNull() const;
    // Check if vector is null

    inline void SetNull();
    // Set all components to 0.

    inline void Set(double xx, double yy, double zz);
    // Assign values to components

    inline void Set(double xx);
    // Assign value to all components

    double Normalize();
    // Normalize to unit this vector

    double Phi() const;
    // The azimuth angle. returns phi from -pi to pi

    double Theta() const;
    // The polar angle.

    inline double CosTheta() const;
    // Cosine of the polar angle.

    inline double Mag2() const;
    // The magnitude squared (rho^2 in spherical coordinate system).

    double Mag() const;
    // The magnitude (rho in spherical coordinate system).

    double Perp2() const;
    // The transverse component (R^2 in cylindrical coordinate system).

    double Perp() const;
    // The transverse component (R in cylindrical coordinate system).

    void RotateX(double);
    // Rotates the vector around the x-axis.

    void RotateY(double);
    // Rotates the vector around the y-axis.

    void RotateZ(double);
    // Rotates the vector around the z-axis.

    inline UVector3& MultiplyByComponents(const UVector3& p);

  public:
    double x;
    double y;
    double z;
};

UVector3 operator + (const UVector3&, const UVector3&);
// Addition of 3-vectors.

UVector3 operator - (const UVector3&, const UVector3&);
// Subtraction of 3-vectors.

double operator * (const UVector3&, const UVector3&);
// Scalar product of 3-vectors.

UVector3 operator * (const UVector3&, double a);
UVector3 operator / (const UVector3&, double a);
UVector3 operator * (double a, const UVector3&);

// Scaling of 3-vectors with a real number

//______________________________________________________________________________
inline UVector3& UVector3::MultiplyByComponents(const UVector3& p)
{
  // Assignment of a UVector3
  x *= p.x;
  y *= p.y;
  z *= p.z;
  return *this;
}

//______________________________________________________________________________
inline UVector3& UVector3::operator = (const UVector3& p)
{
  // Assignment of a UVector3
  if (this == &p)  { return *this; }
  x = p.x;
  y = p.y;
  z = p.z;
  return *this;
}

inline UVector3& UVector3::operator = (const double vect[3])
{
  // Assignment of a C array
  x = vect[0];
  y = vect[1];
  z = vect[2];
  return *this;
}

inline bool UVector3::operator == (const UVector3& v) const
{
  return (v.x == x && v.y == y && v.z == z) ? true : false;
}

inline bool UVector3::operator != (const UVector3& v) const
{
  return (v.x != x || v.y != y || v.z != z) ? true : false;
}

inline UVector3& UVector3::operator += (const UVector3& p)
{
  x += p.x;
  y += p.y;
  z += p.z;
  return *this;
}

inline UVector3& UVector3::operator -= (const UVector3& p)
{
  x -= p.x;
  y -= p.y;
  z -= p.z;
  return *this;
}

inline UVector3 UVector3::operator - () const
{
  return UVector3(-x, -y, -z);
}

inline UVector3& UVector3::operator *= (double a)
{
  x *= a;
  y *= a;
  z *= a;
  return *this;
}

inline UVector3& UVector3::operator /= (double a)
{
  a = 1. / a;
  x *= a;
  y *= a;
  z *= a;
  return *this;
}

inline bool UVector3::IsNull() const
{
  return ((std::abs(x) + std::abs(y) + std::abs(z)) == 0.0) ? true : false;
}

/*
inline void UVector3::SetNull() {
x = y = z = 0.0;
}
*/

inline void UVector3::Set(double xx, double yy, double zz)
{
  x = xx;
  y = yy;
  z = zz;
}

inline void UVector3::Set(double xx)
{
  x = y = z = xx;
}

inline double UVector3::Dot(const UVector3& p) const
{
  return x * p.x + y * p.y + z * p.z;
}

inline UVector3 UVector3::Cross(const UVector3& p) const
{
  return UVector3(y * p.z - p.y * z, z * p.x - p.z * x, x * p.y - p.x * y);
}

inline double UVector3::Mag2() const
{
  return x * x + y * y + z * z;
}

inline double UVector3::Perp2() const
{
  return x * x + y * y;
}

inline double UVector3::CosTheta() const
{
  double ptot = Mag();
  return ptot == 0.0 ? 1.0 : z / ptot;
}


inline double& UVector3::operator[](int index)
{
  switch (index)
  {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      return x;
  }
}

inline double UVector3::operator[](int index) const
{
  //  return operator()(index);

  // TODO: test performance of both versions on Linux
  // => first version is slightly faster
  if (true)
  {
    double vec[3] = {x, y, z};
    return vec[index];
  }

  switch (index)
  {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      return 0;
  }
}

inline std::ostream& operator<< (std::ostream& os, const UVector3& v)
{
  return os << "(" << v.x << "," << v.y << "," << v.z << ")";
}

#endif
