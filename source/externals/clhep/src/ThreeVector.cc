// -*- C++ -*-
// $Id:$
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the Hep3Vector class.
//
// See also ThreeVectorR.cc for implementation of Hep3Vector methods which 
// would couple in all the HepRotation methods.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>
#include <iostream>

namespace CLHEP  {

void Hep3Vector::setMag(double ma) {
  double factor = mag();
  if (factor == 0) {
    std::cerr << "Hep3Vector::setMag() - "
              << "zero vector can't be stretched" << std::endl;
  }else{
    factor = ma/factor;
    setX(x()*factor);
    setY(y()*factor);
    setZ(z()*factor);
  }
}

Hep3Vector & Hep3Vector::rotateUz(const Hep3Vector& NewUzVector) {
  // NewUzVector must be normalized !

  double u1 = NewUzVector.x();
  double u2 = NewUzVector.y();
  double u3 = NewUzVector.z();
  double up = u1*u1 + u2*u2;

  if (up>0) {
      up = std::sqrt(up);
      double px = dx,  py = dy,  pz = dz;
      dx = (u1*u3*px - u2*py)/up + u1*pz;
      dy = (u2*u3*px + u1*py)/up + u2*pz;
      dz =    -up*px +             u3*pz;
    }
  else if (u3 < 0.) { dx = -dx; dz = -dz; }      // phi=0  teta=pi
  else {};
  return *this;
}

double Hep3Vector::pseudoRapidity() const {
  double m1 = mag();
  if ( m1==  0   ) return  0.0;   
  if ( m1==  z() ) return  1.0E72;
  if ( m1== -z() ) return -1.0E72;
  return 0.5*std::log( (m1+z())/(m1-z()) );
}

std::ostream & operator<< (std::ostream & os, const Hep3Vector & v) {
  return os << "(" << v.x() << "," << v.y() << "," << v.z() << ")";
}

void ZMinput3doubles ( std::istream & is, const char * type,
                       double & x, double & y, double & z );

std::istream & operator>>(std::istream & is, Hep3Vector & v) {
  double x, y, z;
  ZMinput3doubles ( is, "Hep3Vector", x, y, z );
  v.set(x, y, z);
  return  is;
}  // operator>>()

const Hep3Vector HepXHat(1.0, 0.0, 0.0);
const Hep3Vector HepYHat(0.0, 1.0, 0.0);
const Hep3Vector HepZHat(0.0, 0.0, 1.0);

//-------------------
//
// New methods introduced when ZOOM PhysicsVectors was merged in:
//
//-------------------

Hep3Vector & Hep3Vector::rotateX (double phi1) {
  double sinphi = std::sin(phi1);
  double cosphi = std::cos(phi1);
  double ty;
  ty = dy * cosphi - dz * sinphi;
  dz = dz * cosphi + dy * sinphi;
  dy = ty;
  return *this;
} /* rotateX */

Hep3Vector & Hep3Vector::rotateY (double phi1) {
  double sinphi = std::sin(phi1);
  double cosphi = std::cos(phi1);
  double tz;
  tz = dz * cosphi - dx * sinphi;
  dx = dx * cosphi + dz * sinphi;
  dz = tz;
  return *this;
} /* rotateY */

Hep3Vector & Hep3Vector::rotateZ (double phi1) {
  double sinphi = std::sin(phi1);
  double cosphi = std::cos(phi1);
  double tx;
  tx = dx * cosphi - dy * sinphi;
  dy = dy * cosphi + dx * sinphi;
  dx = tx;
  return *this;
} /* rotateZ */

bool Hep3Vector::isNear(const Hep3Vector & v, double epsilon) const {
  double limit = dot(v)*epsilon*epsilon;
  return ( (*this - v).mag2() <= limit );
} /* isNear() */

double Hep3Vector::howNear(const Hep3Vector & v ) const {
  // | V1 - V2 | **2  / V1 dot V2, up to 1
  double d   = (*this - v).mag2();
  double vdv = dot(v);
  if ( (vdv > 0) && (d < vdv)  ) {
    return std::sqrt (d/vdv);
  } else if ( (vdv == 0) && (d == 0) ) {
    return 0;
  } else {
    return 1;
  }
} /* howNear */

double Hep3Vector::deltaPhi  (const Hep3Vector & v2) const {
  double dphi = v2.getPhi() - getPhi();
  if ( dphi > CLHEP::pi ) {
    dphi -= CLHEP::twopi;
  } else if ( dphi <= -CLHEP::pi ) {
    dphi += CLHEP::twopi;
  }
  return dphi;
} /* deltaPhi */

double Hep3Vector::deltaR ( const Hep3Vector & v ) const {
  double a = eta() - v.eta();
  double b = deltaPhi(v); 
  return std::sqrt ( a*a + b*b );
} /* deltaR */

double Hep3Vector::cosTheta(const Hep3Vector & q) const {
  double arg;
  double ptot2 = mag2()*q.mag2();
  if(ptot2 <= 0) {
    arg = 0.0;
  }else{
    arg = dot(q)/std::sqrt(ptot2);
    if(arg >  1.0) arg =  1.0;
    if(arg < -1.0) arg = -1.0;
  }
  return arg;
}

double Hep3Vector::cos2Theta(const Hep3Vector & q) const {
  double arg;
  double ptot2 = mag2();
  double qtot2 = q.mag2();
  if ( ptot2 == 0 || qtot2 == 0 )  {
    arg = 1.0;
  }else{
    double pdq = dot(q);
    arg = (pdq/ptot2) * (pdq/qtot2);
        // More naive methods overflow on vectors which can be squared
        // but can't be raised to the 4th power.
    if(arg >  1.0) arg =  1.0;
 }
 return arg;
}

void Hep3Vector::setEta (double eta1) {
  double phi1 = 0;
  double r1;
  if ( (dx == 0) && (dy == 0) ) {
    if (dz == 0) {
      std::cerr << "Hep3Vector::setEta() - "
                << "Attempt to set eta of zero vector -- vector is unchanged"
                << std::endl;
      return;
    }
  std::cerr << "Hep3Vector::setEta() - "
            << "Attempt to set eta of vector along Z axis -- will use phi = 0"
            << std::endl;
    r1 = std::fabs(dz);
  } else {
    r1 = getR();
    phi1 = getPhi();
  }
  double tanHalfTheta = std::exp ( -eta1 );
  double cosTheta1 =
        (1 - tanHalfTheta*tanHalfTheta) / (1 + tanHalfTheta*tanHalfTheta);
  dz = r1 * cosTheta1;
  double rho1 = r1*std::sqrt(1 - cosTheta1*cosTheta1);
  dy = rho1 * std::sin (phi1);
  dx = rho1 * std::cos (phi1);
  return;
}

void Hep3Vector::setCylTheta (double theta1) {

  // In cylindrical coords, set theta while keeping rho and phi fixed

  if ( (dx == 0) && (dy == 0) ) {
    if (dz == 0) {
      std::cerr << "Hep3Vector::setCylTheta() - "
                << "Attempt to set cylTheta of zero vector -- vector is unchanged"
                << std::endl;
      return;
    }
    if (theta1 == 0) {
      dz = std::fabs(dz);
      return;
    }
    if (theta1 == CLHEP::pi) {
      dz = -std::fabs(dz);
      return;
    }
    std::cerr << "Hep3Vector::setCylTheta() - "
      << "Attempt set cylindrical theta of vector along Z axis "
      << "to a non-trivial value, while keeping rho fixed -- "
      << "will return zero vector" << std::endl;
    dz = 0;
    return;
  }
  if ( (theta1 < 0) || (theta1 > CLHEP::pi) ) {
    std::cerr << "Hep3Vector::setCylTheta() - "
      << "Setting Cyl theta of a vector based on a value not in [0, PI]"
      << std::endl;
        // No special return needed if warning is ignored.
  }
  double phi1 (getPhi());
  double rho1 = getRho();
  if ( (theta1 == 0) || (theta1 == CLHEP::pi) ) {
    std::cerr << "Hep3Vector::setCylTheta() - "
      << "Attempt to set cylindrical theta to 0 or PI "
      << "while keeping rho fixed -- infinite Z will be computed"
      << std::endl;
      dz = (theta1==0) ? 1.0E72 : -1.0E72;
    return;
  }
  dz = rho1 / std::tan (theta1);
  dy = rho1 * std::sin (phi1);
  dx = rho1 * std::cos (phi1);

} /* setCylTheta */

void Hep3Vector::setCylEta (double eta1) {

  // In cylindrical coords, set eta while keeping rho and phi fixed

  double theta1 = 2 * std::atan ( std::exp (-eta1) );

        //-| The remaining code is similar to setCylTheta,  The reason for
        //-| using a copy is so as to be able to change the messages in the
        //-| ZMthrows to say eta rather than theta.  Besides, we assumedly
        //-| need not check for theta of 0 or PI.

  if ( (dx == 0) && (dy == 0) ) {
    if (dz == 0) {
      std::cerr << "Hep3Vector::setCylEta() - "
        << "Attempt to set cylEta of zero vector -- vector is unchanged"
        << std::endl;
      return;
    }
    if (theta1 == 0) {
      dz = std::fabs(dz);
      return;
    }
    if (theta1 == CLHEP::pi) {
      dz = -std::fabs(dz);
      return;
    }
    std::cerr << "Hep3Vector::setCylEta() - "
      << "Attempt set cylindrical eta of vector along Z axis "
      << "to a non-trivial value, while keeping rho fixed -- "
      << "will return zero vector" << std::endl;
    dz = 0;
    return;
  }
  double phi1 (getPhi());
  double rho1 = getRho();
  dz = rho1 / std::tan (theta1);
  dy = rho1 * std::sin (phi1);
  dx = rho1 * std::cos (phi1);

} /* setCylEta */


Hep3Vector operator/  ( const Hep3Vector & v1, double c ) {
//  if (c == 0) {
//    std::cerr << "Hep3Vector::operator/ () - "
//      << "Attempt to divide vector by 0 -- "
//      << "will produce infinities and/or NANs" << std::endl;
//  } 
  double   oneOverC = 1.0/c;
  return Hep3Vector  (  v1.x() * oneOverC,
                        v1.y() * oneOverC,
                        v1.z() * oneOverC );
} /* v / c */

Hep3Vector & Hep3Vector::operator/= (double c) {
//  if (c == 0) {
//    std::cerr << "Hep3Vector::operator/ () - "
//      << "Attempt to do vector /= 0 -- "
//      << "division by zero would produce infinite or NAN components"
//      << std::endl;
//  }
  double oneOverC = 1.0/c;
  dx *= oneOverC;
  dy *= oneOverC;
  dz *= oneOverC;
  return *this;
}

double Hep3Vector::tolerance = Hep3Vector::ToleranceTicks * 2.22045e-16;

}  // namespace CLHEP
