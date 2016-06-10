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
// UVector2
//
// Class description:
//
//    UVector2 is a general 2-vector class defining vectors in two
//    dimension using double components.
//
// 19.09.12 Marek Gayer
//          Created from original implementation in CLHEP
// --------------------------------------------------------------------

#ifndef UVECTOR2_H
#define UVECTOR2_H

#include <cmath>
#include <iostream>

#include "UVector3.hh"

// Declarations of classes and global methods
class UVector2;
std::ostream& operator << (std::ostream&, const UVector2&);
//std::istream & operator >> (std::istream &, UVector2 &);
inline double operator * (const UVector2& a, const UVector2& b);
inline UVector2 operator * (const UVector2& p, double a);
inline UVector2 operator * (double a, const UVector2& p);
UVector2 operator / (const UVector2& p, double a);
inline UVector2 operator + (const UVector2& a, const UVector2& b);
inline UVector2 operator - (const UVector2& a, const UVector2& b);

/**
* @author
* @ingroup vector
*/
class UVector2
{

  public:

    enum { X = 0, Y = 1, NUM_COORDINATES = 2, SIZE = NUM_COORDINATES };
    // Safe indexing of the coordinates when using with matrices, arrays, etc.

    inline UVector2(double x = 0.0, double y = 0.0);
    // The constructor.

    inline UVector2(const UVector2& p);
    // The copy constructor.

    explicit UVector2(const UVector3& s);
    // "demotion" constructor"
    // WARNING -- THIS IGNORES THE Z COMPONENT OF THE UVector3.
    //    SO IN GENERAL, UVector2(v)==v WILL NOT HOLD!

    inline ~UVector2();
    // The destructor.

//  inline double x() const;
//  inline double y() const;
    // The components in cartesian coordinate system.

    double operator()(int i) const;
    inline double operator [](int i) const;
    // Get components by index.  0-based.

    double& operator()(int i);
    inline double& operator [](int i);
    // Set components by index.  0-based.

    inline void setX(double x);
    inline void setY(double y);
    inline void set(double x, double y);
    // Set the components in cartesian coordinate system.

    inline double phi() const;
    // The azimuth angle.

    inline double mag2() const;
    // The magnitude squared.

    inline double mag() const;
    // The magnitude.

    inline double r() const;
    // r in polar coordinates (r, phi):  equal to mag().

    inline void setPhi(double phi);
    // Set phi keeping mag constant.

    inline void setMag(double r);
    // Set magnitude keeping phi constant.

    inline void setR(double r);
    // Set R keeping phi constant.  Same as setMag.

    inline void setPolar(double r, double phi);
    // Set by polar coordinates.

    inline UVector2& operator = (const UVector2& p);
    // Assignment.

    inline bool operator == (const UVector2& v) const;
    inline bool operator != (const UVector2& v) const;
    // Comparisons.

    int compare(const UVector2& v) const;
    bool operator > (const UVector2& v) const;
    bool operator < (const UVector2& v) const;
    bool operator>= (const UVector2& v) const;
    bool operator<= (const UVector2& v) const;
    // dictionary ordering according to y, then x component

    static inline double getTolerance();
    static double setTolerance(double tol);

    double howNear(const UVector2& p) const;
    bool isNear(const UVector2& p, double epsilon = tolerance) const;

    double howParallel(const UVector2& p) const;
    bool isParallel
    (const UVector2& p, double epsilon = tolerance) const;

    double howOrthogonal(const UVector2& p) const;
    bool isOrthogonal
    (const UVector2& p, double epsilon = tolerance) const;

    inline UVector2& operator += (const UVector2& p);
    // Addition.

    inline UVector2& operator -= (const UVector2& p);
    // Subtraction.

    inline UVector2 operator - () const;
    // Unary minus.

    inline UVector2& operator *= (double a);
    // Scaling with real numbers.

    inline UVector2 unit() const;
    // Unit vector parallel to this.

    inline UVector2 orthogonal() const;
    // Vector orthogonal to this.

    inline double dot(const UVector2& p) const;
    // Scalar product.

    inline double angle(const UVector2&) const;
    // The angle w.r.t. another 2-vector.

    void rotate(double);
    // Rotates the UVector2.

    operator UVector3() const;
    // Cast a UVector2 as a UVector3.

    // The remaining methods are friends, thus defined at global scope:
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    friend std::ostream& operator<< (std::ostream&, const UVector2&);
    // Output to a stream.

    inline friend double operator * (const UVector2& a,
                                     const UVector2& b);
    // Scalar product.

    inline friend UVector2 operator * (const UVector2& p, double a);
    // v*c

    inline friend UVector2 operator * (double a, const UVector2& p);
    // c*v

    friend UVector2 operator / (const UVector2& p, double a);
    // v/c

    inline friend UVector2 operator + (const UVector2& a,
                                       const UVector2& b);
    // v1+v2

    inline friend UVector2 operator - (const UVector2& a,
                                       const UVector2& b);
    // v1-v2

    enum { ZMpvToleranceTicks = 100 };

    double x;
    double y;
    // The components.

  private:

    static double tolerance;
    // default tolerance criterion for isNear() to return true.

};  // UVector2

static const UVector2 X_HAT2(1.0, 0.0);
static const UVector2 Y_HAT2(0.0, 1.0);


/*
inline double UVector2::x() const {
  return x;
}

inline double UVector2::y() const {
  return y;
}
*/

inline UVector2::UVector2(double x1, double y1)
  : x(x1), y(y1) {}

inline UVector2::UVector2(const UVector3& s1)
  : x(s1.x()), y(s1.y()) {}

inline void UVector2::setX(double x1)
{
  x = x1;
}

inline void UVector2::setY(double y1)
{
  y = y1;
}

inline void UVector2::set(double x1, double y1)
{
  x = x1;
  y = y1;
}

double& UVector2::operator[](int i)
{
  return operator()(i);
}
double   UVector2::operator[](int i) const
{
  return operator()(i);
}

inline UVector2::UVector2(const UVector2& p)
  : x(p.x), y(p.y) {}

inline UVector2::~UVector2() {}

inline UVector2& UVector2::operator = (const UVector2& p)
{
  if (this == &p)  { return *this; }
  x = p.x;
  y = p.y;
  return *this;
}

inline bool UVector2::operator == (const UVector2& v) const
{
  return (v.x == x && v.y == y) ? true : false;
}

inline bool UVector2::operator != (const UVector2& v) const
{
  return (v.x != x || v.y != y) ? true : false;
}

inline UVector2& UVector2::operator += (const UVector2& p)
{
  x += p.x;
  y += p.y;
  return *this;
}

inline UVector2& UVector2::operator -= (const UVector2& p)
{
  x -= p.x;
  y -= p.y;
  return *this;
}

inline UVector2 UVector2::operator - () const
{
  return UVector2(-x, -y);
}

inline UVector2& UVector2::operator *= (double a)
{
  x *= a;
  y *= a;
  return *this;
}

inline double UVector2::dot(const UVector2& p) const
{
  return x * p.x + y * p.y;
}

inline double UVector2::mag2() const
{
  return x * x + y * y;
}

inline double UVector2::mag() const
{
  return std::sqrt(mag2());
}

inline double UVector2::r() const
{
  return std::sqrt(mag2());
}

inline UVector2 UVector2::unit() const
{
  double tot = mag2();
  UVector2 p(*this);
  return tot > 0.0 ? p *= (1.0 / std::sqrt(tot)) : UVector2(1, 0);
}

inline UVector2 UVector2::orthogonal() const
{
  double x1 = std::fabs(x), y1 = std::fabs(y);
  if (x1 < y1)
  {
    return UVector2(y, -x);
  }
  else
  {
    return UVector2(-y, x);
  }
}

inline double UVector2::phi() const
{
  return x == 0.0 && y == 0.0 ? 0.0 : std::atan2(y, x);
}

inline double UVector2::angle(const UVector2& q) const
{
  double ptot2 = mag2() * q.mag2();
  return ptot2 <= 0.0 ? 0.0 : std::acos(dot(q) / std::sqrt(ptot2));
}

inline void UVector2::setMag(double r1)
{
  double ph = phi();
  setX(r1 * std::cos(ph));
  setY(r1 * std::sin(ph));
}

inline void UVector2::setR(double r1)
{
  setMag(r1);
}

inline void UVector2::setPhi(double phi1)
{
  double ma = mag();
  setX(ma * std::cos(phi1));
  setY(ma * std::sin(phi1));
}

inline void UVector2::setPolar(double r1, double phi1)
{
  setX(r1 * std::cos(phi1));
  setY(r1 * std::sin(phi1));
}

inline UVector2 operator + (const UVector2& a, const UVector2& b)
{
  return UVector2(a.x + b.x, a.y + b.y);
}

inline UVector2 operator - (const UVector2& a, const UVector2& b)
{
  return UVector2(a.x - b.x, a.y - b.y);
}

inline UVector2 operator * (const UVector2& p, double a)
{
  return UVector2(a * p.x, a * p.y);
}

inline UVector2 operator * (double a, const UVector2& p)
{
  return UVector2(a * p.x, a * p.y);
}

inline double operator * (const UVector2& a, const UVector2& b)
{
  return a.dot(b);
}

inline double UVector2::getTolerance()
{
  return tolerance;
}


#endif /* UVECTOR2_H */
