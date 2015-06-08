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
      x_ = y_ = z_ = 0.0;
    }
    UVector3(double xval, double yval, double zval)
    {
      x_ = xval;
      y_ = yval;
      z_ = zval;
    }
    UVector3(double theta, double phi);
    UVector3(const double coord[3])
    {
      x_ = coord[0];
      y_ = coord[1];
      z_ = coord[2];
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

    inline double& x();

    inline double x() const;

    inline double& y();

    inline double y() const;

    inline double& z();

    inline double z() const;

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

  private:
    double x_;
    double y_;
    double z_;
};


inline UVector3 operator + (const UVector3& a, const UVector3& b)
{
  return UVector3(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
}

inline UVector3 operator - (const UVector3& a, const UVector3& b)
{
  return UVector3(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
}

inline UVector3 operator * (const UVector3& p, double a)
{
  return UVector3(a * p.x(), a * p.y(), a * p.z());
}

inline UVector3 operator / (const UVector3& p, double a)
{
  a = 1. / a;
  return UVector3(a * p.x(), a * p.y(), a * p.z());
}

inline UVector3 operator * (double a, const UVector3& p)
{
  return UVector3(a * p.x(), a * p.y(), a * p.z());
}


//______________________________________________________________________________
inline UVector3& UVector3::MultiplyByComponents(const UVector3& p)
{
  // Assignment of a UVector3
  x_ *= p.x_;
  y_ *= p.y_;
  z_ *= p.z_;
  return *this;
}

//______________________________________________________________________________
inline UVector3& UVector3::operator = (const UVector3& p)
{
  // Assignment of a UVector3
  if (this == &p)  { return *this; }
  x_ = p.x_;
  y_ = p.y_;
  z_ = p.z_;
  return *this;
}

inline UVector3& UVector3::operator = (const double vect[3])
{
  // Assignment of a C array
  x_ = vect[0];
  y_ = vect[1];
  z_ = vect[2];
  return *this;
}

inline bool UVector3::operator == (const UVector3& v) const
{
  return (v.x_ == x_ && v.y_ == y_ && v.z_ == z_) ? true : false;
}

inline bool UVector3::operator != (const UVector3& v) const
{
  return (v.x_ != x_ || v.y_ != y_ || v.z_ != z_) ? true : false;
}

inline UVector3& UVector3::operator += (const UVector3& p)
{
  x_ += p.x_;
  y_ += p.y_;
  z_ += p.z_;
  return *this;
}

inline UVector3& UVector3::operator -= (const UVector3& p)
{
  x_ -= p.x_;
  y_ -= p.y_;
  z_ -= p.z_;
  return *this;
}

inline UVector3 UVector3::operator - () const
{
  return UVector3(-x_, -y_, -z_);
}

inline UVector3& UVector3::operator *= (double a)
{
  x_ *= a;
  y_ *= a;
  z_ *= a;
  return *this;
}

inline UVector3& UVector3::operator /= (double a)
{
  a = 1. / a;
  x_ *= a;
  y_ *= a;
  z_ *= a;
  return *this;
}

inline bool UVector3::IsNull() const
{
  return ((std::abs(x_) + std::abs(y_) + std::abs(z_)) == 0.0) ? true : false;
}

inline void UVector3::Set(double xx, double yy, double zz)
{
  x_ = xx;
  y_ = yy;
  z_ = zz;
}

inline void UVector3::Set(double xx)
{
  x_ = y_ = z_ = xx;
}

inline double UVector3::Dot(const UVector3& p) const
{
  return x_ * p.x_ + y_ * p.y_ + z_ * p.z_;
}

inline UVector3 UVector3::Cross(const UVector3& p) const
{
  return UVector3(y_ * p.z_ - p.y_ * z_, z_ * p.x_ - p.z_ * x_,
                  x_ * p.y_ - p.x_ * y_);
}

inline double UVector3::Mag2() const
{
  return x_ * x_ + y_ * y_ + z_ * z_;
}

inline double UVector3::Perp2() const
{
  return x_ * x_ + y_ * y_;
}

inline double UVector3::CosTheta() const
{
  double ptot = Mag();
  return ptot == 0.0 ? 1.0 : z_ / ptot;
}


inline double& UVector3::operator[](int index)
{
  switch (index)
  {
    case 0:
      return x_;
    case 1:
      return y_;
    case 2:
      return z_;
    default:
      return x_;
  }
}

inline double UVector3::operator[](int index) const
{
  //  return operator()(index);

  // TODO: test performance of both versions on Linux
  // => first version is slightly faster
  if (true)
  {
    double vec[3] = {x_, y_, z_};
    return vec[index];
  }

  switch (index)
  {
    case 0:
      return x_;
    case 1:
      return y_;
    case 2:
      return z_;
    default:
      return 0;
  }
}

inline double& UVector3::x() { return x_; }

inline double UVector3::x() const { return x_; }

inline double& UVector3::y() { return y_; }

inline double UVector3::y() const { return y_; }

inline double& UVector3::z() { return z_; }

inline double UVector3::z() const { return z_; }

inline std::ostream& operator<< (std::ostream& os, const UVector3& v)
{
  return os << "(" << v.x() << "," << v.y() << "," << v.z() << ")";
}

#endif
