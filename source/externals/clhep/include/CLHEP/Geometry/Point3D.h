// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// History:
// 09.09.96 E.Chernyaev - initial version
// 12.06.01 E.Chernyaev - CLHEP-1.7: introduction of BasicVector3D to decouple
//                        the functionality from CLHEP::Hep3Vector
// 01.04.03 E.Chernyaev - CLHEP-1.9: template version
//

#ifndef HEP_POINT3D_H
#define HEP_POINT3D_H

#include <iosfwd>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/BasicVector3D.h"

namespace HepGeom {

  class Transform3D;

  /**
   * Geometrical 3D Point.
   * This is just a declaration of the class needed to define
   * specializations Point3D<float> and Point3D<double>.
   *
   * @ingroup geometry
   * @author Evgeni Chernyaev <Evgueni.Tcherniaev@cern.ch>
   */
  template<class T>
  class Point3D : public BasicVector3D<T> {};

  /**
   * Geometrical 3D Point with components of float type.
   *
   * @author Evgeni Chernyaev <Evgueni.Tcherniaev@cern.ch>
   * @ingroup geometry
   */
  template<>
  class Point3D<float> : public BasicVector3D<float> {
  public:
    /**
     * Default constructor. */
    Point3D() = default;

    /**
     * Constructor from three numbers. */
    Point3D(float x1, float y1, float z1) : BasicVector3D<float>(x1,y1,z1) {}

    /**
     * Constructor from array of floats. */
    explicit Point3D(const float * a)
      : BasicVector3D<float>(a[0],a[1],a[2]) {}

    /**
     * Copy constructor. */
    Point3D(const Point3D<float> &) = default;

    /**
     * Move constructor. */
    Point3D(Point3D<float> &&) = default;

    /**
     * Constructor from BasicVector3D<float>. */
    Point3D(const BasicVector3D<float> & v) : BasicVector3D<float>(v) {}

    /**
     * Destructor. */
    ~Point3D() = default;

    /**
     * Assignment. */
    Point3D<float> & operator=(const Point3D<float> &) = default;

    /**
     * Assignment from BasicVector3D<float>. */
    Point3D<float> & operator=(const BasicVector3D<float> & v) {
      this->BasicVector3D<float>::operator=(v);
      return *this;
    }

    /**
     * Move assignment. */
    Point3D<float> & operator=(Point3D<float> &&) = default;

    /**
     * Returns distance to the origin squared. */
    float distance2() const { return mag2(); }

    /**
     * Returns distance to the point squared. */
    float distance2(const Point3D<float> & p) const {
      float dx = p.x()-x(), dy = p.y()-y(), dz = p.z()-z();
      return dx*dx + dy*dy + dz*dz;
    }

    /**
     * Returns distance to the origin. */
    float distance() const { return std::sqrt(distance2()); }

    /**
     * Returns distance to the point. */
    float distance(const Point3D<float> & p) const {
      return std::sqrt(distance2(p));
    }

    /**
     * Transformation by Transform3D. */
    Point3D<float> & transform(const Transform3D & m);
  };

  /**
   * Transformation of Point3D<float> by Transform3D.
   * @relates Point3D
   */
  Point3D<float>
  operator*(const Transform3D & m, const Point3D<float> & p);

  /**
   * Geometrical 3D Point with components of double type.
   *
   * @author Evgeni Chernyaev <Evgueni.Tcherniaev@cern.ch>
   * @ingroup geometry
   */
  template<>
  class Point3D<double> : public BasicVector3D<double> {
  public:
    /**
     * Default constructor. */
    Point3D() = default;

    /**
     * Constructor from three numbers. */
    Point3D(double x1, double y1, double z1) : BasicVector3D<double>(x1,y1,z1) {}

    /**
     * Constructor from array of floats. */
    explicit Point3D(const float * a)
      : BasicVector3D<double>(a[0],a[1],a[2]) {}

    /**
     * Constructor from array of doubles. */
    explicit Point3D(const double * a)
      : BasicVector3D<double>(a[0],a[1],a[2]) {}

    /**
     * Copy constructor. */
    Point3D(const Point3D<double> &) = default;

    /**
     * Move constructor. */
    Point3D(Point3D<double> &&) = default;

    /**
     * Constructor from BasicVector3D<float>. */
    Point3D(const BasicVector3D<float> & v) : BasicVector3D<double>(v) {}

    /**
     * Constructor from BasicVector3D<double>. */
    Point3D(const BasicVector3D<double> & v) : BasicVector3D<double>(v) {}

    /**
     * Destructor. */
    ~Point3D() = default;

    /**
     * Constructor from CLHEP::Hep3Vector.
     * This constructor is needed only for backward compatibility and
     * in principle should be absent.
     */
    Point3D(const CLHEP::Hep3Vector & v)
      : BasicVector3D<double>(v.x(),v.y(),v.z()) {}

    /**
     * Conversion (cast) to CLHEP::Hep3Vector.
     * This operator is needed only for backward compatibility and
     * in principle should not exit.
     */
    operator CLHEP::Hep3Vector () const { return CLHEP::Hep3Vector(x(),y(),z()); }

    /**
     * Assignment. */
    Point3D<double> & operator=(const Point3D<double> &) = default;

    /**
     * Assignment from BasicVector3D<float>. */
    Point3D<double> & operator=(const BasicVector3D<float> & v) {
      this->BasicVector3D<double>::operator=(v);
      return *this;
    }

    /**
     * Assignment from BasicVector3D<double>. */
    Point3D<double> & operator=(const BasicVector3D<double> & v) {
      this->BasicVector3D<double>::operator=(v);
      return *this;
    }

    /**
     * Move assignment. */
    Point3D<double> & operator=(Point3D<double> &&) = default;

    /**
     * Returns distance to the origin squared. */
    double distance2() const { return mag2(); }

    /**
     * Returns distance to the point squared. */
    double distance2(const Point3D<double> & p) const {
      double dx = p.x()-x(), dy = p.y()-y(), dz = p.z()-z();
      return dx*dx + dy*dy + dz*dz;
    }

    /**
     * Returns distance to the origin. */
    double distance() const { return std::sqrt(distance2()); }

    /**
     * Returns distance to the point. */
    double distance(const Point3D<double> & p) const {
      return std::sqrt(distance2(p));
    }

    /**
     * Transformation by Transform3D. */
    Point3D<double> & transform(const Transform3D & m);
  };

  /**
   * Transformation of Point3D<double> by Transform3D.
   * @relates Point3D
   */
  Point3D<double>
  operator*(const Transform3D & m, const Point3D<double> & p);

} /* namespace HepGeom */

#endif /* HEP_POINT3D_H */
