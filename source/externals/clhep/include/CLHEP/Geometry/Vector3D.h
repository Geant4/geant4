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

#ifndef HEP_VECTOR3D_H
#define HEP_VECTOR3D_H

#include <iosfwd>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/BasicVector3D.h"

namespace HepGeom {

  class Transform3D;

  /**
   * Geometrical 3D Vector.
   * This is just a declaration of the class needed to define
   * specializations Vector3D<float> and Vector3D<double>.
   *
   * @ingroup geometry
   * @author Evgeni Chernyaev <Evgueni.Tcherniaev@cern.ch>
   */
  template<class T>
  class Vector3D : public BasicVector3D<T> {};

  /**
   * Geometrical 3D Vector with components of float type.
   *
   * @author Evgeni Chernyaev <Evgueni.Tcherniaev@cern.ch>
   * @ingroup geometry
   */
  template<>
  class Vector3D<float> : public BasicVector3D<float> {
  public:
    /**
     * Default constructor. */
    Vector3D() = default;

    /**
     * Constructor from three numbers. */
    Vector3D(float x1, float y1, float z1) : BasicVector3D<float>(x1,y1,z1) {}

    /**
     * Constructor from array of floats. */
    explicit Vector3D(const float * a)
      : BasicVector3D<float>(a[0],a[1],a[2]) {}

    /**
     * Copy constructor. */
    Vector3D(const Vector3D<float> &) = default;

    /**
     * Move constructor. */
    Vector3D(Vector3D<float> &&) = default;

    /**
     * Constructor from BasicVector3D<float>. */
    Vector3D(const BasicVector3D<float> & v) : BasicVector3D<float>(v) {}

    /**
     * Destructor. */
    ~Vector3D() = default;

    /**
     * Assignment. */
    Vector3D<float> & operator=(const Vector3D<float> &) = default;

    /**
     * Assignment from BasicVector3D<float>. */
    Vector3D<float> & operator=(const BasicVector3D<float> & v) {
      this->BasicVector3D<float>::operator=(v);
      return *this;
    }

    /**
     * Move assignment. */
    Vector3D<float> & operator=(Vector3D<float> &&) = default;

    /**
     * Transformation by Transform3D. */
    Vector3D<float> & transform(const Transform3D & m);
  };

  /**
   * Transformation of Vector<float> by Transform3D.
   * @relates Vector3D
   */
  Vector3D<float>
  operator*(const Transform3D & m, const Vector3D<float> & v);

  /**
   * Geometrical 3D Vector with components of double type.
   *
   * @author Evgeni Chernyaev <Evgueni.Tcherniaev@cern.ch>
   * @ingroup geometry
   */
  template<>
  class Vector3D<double> : public BasicVector3D<double> {
  public:
    /**
     * Default constructor. */
    Vector3D() = default;

    /**
     * Constructor from three numbers. */
    Vector3D(double x1, double y1, double z1) : BasicVector3D<double>(x1,y1,z1) {}

    /**
     * Constructor from array of floats. */
    explicit Vector3D(const float * a)
      : BasicVector3D<double>(a[0],a[1],a[2]) {}

    /**
     * Constructor from array of doubles. */
    explicit Vector3D(const double * a)
      : BasicVector3D<double>(a[0],a[1],a[2]) {}

    /**
     * Copy constructor. */
    Vector3D(const Vector3D<double> &) = default;

    /**
     * Move constructor. */
    Vector3D(Vector3D<double> &&) = default;

    /**
     * Constructor from BasicVector3D<float>. */
    Vector3D(const BasicVector3D<float> & v) : BasicVector3D<double>(v) {}

    /**
     * Constructor from BasicVector3D<double>. */
    Vector3D(const BasicVector3D<double> & v) : BasicVector3D<double>(v) {}

    /**
     * Destructor. */
    ~Vector3D() = default;

    /**
     * Constructor from CLHEP::Hep3Vector.
     * This constructor is needed only for backward compatibility and
     * in principle should be absent.
     */
    Vector3D(const CLHEP::Hep3Vector & v)
      : BasicVector3D<double>(v.x(),v.y(),v.z()) {}

    /**
     * Conversion (cast) to CLHEP::Hep3Vector.
     * This operator is needed only for backward compatibility and
     * in principle should not exit.
     */
    operator CLHEP::Hep3Vector () const { return CLHEP::Hep3Vector(x(),y(),z()); }

    /**
     * Assignment. */
    Vector3D<double> & operator=(const Vector3D<double> &) = default;

    /**
     * Assignment from BasicVector3D<float>. */
    Vector3D<double> & operator=(const BasicVector3D<float> & v) {
      this->BasicVector3D<double>::operator=(v);
      return *this;
    }

    /**
     * Assignment from BasicVector3D<double>. */
    Vector3D<double> & operator=(const BasicVector3D<double> & v) {
      this->BasicVector3D<double>::operator=(v);
      return *this;
    }

    /**
     * Move assignment. */
    Vector3D<double> & operator=(Vector3D<double> &&) = default;

    /**
     * Transformation by Transform3D. */
    Vector3D<double> & transform(const Transform3D & m);
  };

  /**
   * Transformation of Vector<double> by Transform3D.
   * @relates Vector3D
   */
  Vector3D<double>
  operator*(const Transform3D & m, const Vector3D<double> & v);

} /* namespace HepGeom */

#endif /* HEP_VECTOR3D_H */
