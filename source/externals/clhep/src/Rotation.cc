// -*- C++ -*-
// $Id:$
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the parts of the the HepRotation class which
// were present in the original CLHEP before the merge with ZOOM PhysicsVectors.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <iostream>
#include <cmath>

namespace CLHEP  {

static inline double safe_acos (double x) {
  if (std::abs(x) <= 1.0) return std::acos(x);
  return ( (x>0) ? 0 : CLHEP::pi );
}

double HepRotation::operator() (int i, int j) const {
  if (i == 0) {
    if (j == 0) { return xx(); }
    if (j == 1) { return xy(); }
    if (j == 2) { return xz(); } 
  } else if (i == 1) {
    if (j == 0) { return yx(); }
    if (j == 1) { return yy(); }
    if (j == 2) { return yz(); } 
  } else if (i == 2) {
    if (j == 0) { return zx(); }
    if (j == 1) { return zy(); }
    if (j == 2) { return zz(); } 
  } 
  std::cerr << "HepRotation subscripting: bad indices "
       << "(" << i << "," << j << ")" << std::endl;
  return 0.0;
} 

HepRotation & HepRotation::rotate(double a, const Hep3Vector& aaxis) {
  if (a != 0.0) {
    double ll = aaxis.mag();
    if (ll == 0.0) {
      std::cerr << "HepRotation::rotate() - "
                << "HepRotation: zero axis" << std::endl;
    }else{
      double sa = std::sin(a), ca = std::cos(a);
      double dx = aaxis.x()/ll, dy = aaxis.y()/ll, dz = aaxis.z()/ll;   
      HepRotation m1(
	ca+(1-ca)*dx*dx,          (1-ca)*dx*dy-sa*dz,    (1-ca)*dx*dz+sa*dy,
	   (1-ca)*dy*dx+sa*dz, ca+(1-ca)*dy*dy,          (1-ca)*dy*dz-sa*dx,
	   (1-ca)*dz*dx-sa*dy,    (1-ca)*dz*dy+sa*dx, ca+(1-ca)*dz*dz );
      transform(m1);
    }
  }
  return *this;
}

HepRotation & HepRotation::rotateX(double a) {
  double c1 = std::cos(a);
  double s1 = std::sin(a);
  double x1 = ryx, y1 = ryy, z1 = ryz; 
  ryx = c1*x1 - s1*rzx;
  ryy = c1*y1 - s1*rzy;
  ryz = c1*z1 - s1*rzz;
  rzx = s1*x1 + c1*rzx;
  rzy = s1*y1 + c1*rzy;
  rzz = s1*z1 + c1*rzz;
  return *this;
}

HepRotation & HepRotation::rotateY(double a){
  double c1 = std::cos(a);
  double s1 = std::sin(a);
  double x1 = rzx, y1 = rzy, z1 = rzz; 
  rzx = c1*x1 - s1*rxx;
  rzy = c1*y1 - s1*rxy;
  rzz = c1*z1 - s1*rxz;
  rxx = s1*x1 + c1*rxx;
  rxy = s1*y1 + c1*rxy;
  rxz = s1*z1 + c1*rxz;
  return *this;
}

HepRotation & HepRotation::rotateZ(double a) {
  double c1 = std::cos(a);
  double s1 = std::sin(a);
  double x1 = rxx, y1 = rxy, z1 = rxz; 
  rxx = c1*x1 - s1*ryx;
  rxy = c1*y1 - s1*ryy;
  rxz = c1*z1 - s1*ryz;
  ryx = s1*x1 + c1*ryx;
  ryy = s1*y1 + c1*ryy;
  ryz = s1*z1 + c1*ryz;
  return *this;
}

HepRotation & HepRotation::rotateAxes(const Hep3Vector &newX,
				      const Hep3Vector &newY,
				      const Hep3Vector &newZ) {
  double del = 0.001;
  Hep3Vector w = newX.cross(newY);

  if (std::abs(newZ.x()-w.x()) > del ||
      std::abs(newZ.y()-w.y()) > del ||
      std::abs(newZ.z()-w.z()) > del ||
      std::abs(newX.mag2()-1.) > del ||
      std::abs(newY.mag2()-1.) > del || 
      std::abs(newZ.mag2()-1.) > del ||
      std::abs(newX.dot(newY)) > del ||
      std::abs(newY.dot(newZ)) > del ||
      std::abs(newZ.dot(newX)) > del) {
    std::cerr << "HepRotation::rotateAxes: bad axis vectors" << std::endl;
    return *this;
  }else{
    return transform(HepRotation(newX.x(), newY.x(), newZ.x(),
                                 newX.y(), newY.y(), newZ.y(),
                                 newX.z(), newY.z(), newZ.z()));
  }
}

double HepRotation::phiX() const {
  return (yx() == 0.0 && xx() == 0.0) ? 0.0 : std::atan2(yx(),xx());
}

double HepRotation::phiY() const {
  return (yy() == 0.0 && xy() == 0.0) ? 0.0 : std::atan2(yy(),xy());
}

double HepRotation::phiZ() const {
  return (yz() == 0.0 && xz() == 0.0) ? 0.0 : std::atan2(yz(),xz());
}

double HepRotation::thetaX() const {
  return safe_acos(zx());
}

double HepRotation::thetaY() const {
  return safe_acos(zy());
}

double HepRotation::thetaZ() const {
  return safe_acos(zz());
}

void HepRotation::getAngleAxis(double &angle, Hep3Vector &aaxis) const {
  double cosa  = 0.5*(xx()+yy()+zz()-1);
  double cosa1 = 1-cosa;
  if (cosa1 <= 0) {
    angle = 0;
    aaxis  = Hep3Vector(0,0,1);
  }else{
    double x=0, y=0, z=0;
    if (xx() > cosa) x = std::sqrt((xx()-cosa)/cosa1);
    if (yy() > cosa) y = std::sqrt((yy()-cosa)/cosa1);
    if (zz() > cosa) z = std::sqrt((zz()-cosa)/cosa1);
    if (zy() < yz()) x = -x;
    if (xz() < zx()) y = -y;
    if (yx() < xy()) z = -z;
    angle = (cosa < -1.) ? std::acos(-1.) : std::acos(cosa);
    aaxis  = Hep3Vector(x,y,z);
  }
}

bool HepRotation::isIdentity() const {
  return  (rxx == 1.0 && rxy == 0.0 && rxz == 0.0 &&
           ryx == 0.0 && ryy == 1.0 && ryz == 0.0 &&
           rzx == 0.0 && rzy == 0.0 && rzz == 1.0) ? true : false;
}

int HepRotation::compare ( const HepRotation & r ) const {
       if (rzz<r.rzz) return -1; else if (rzz>r.rzz) return 1;
  else if (rzy<r.rzy) return -1; else if (rzy>r.rzy) return 1;
  else if (rzx<r.rzx) return -1; else if (rzx>r.rzx) return 1;
  else if (ryz<r.ryz) return -1; else if (ryz>r.ryz) return 1;
  else if (ryy<r.ryy) return -1; else if (ryy>r.ryy) return 1;
  else if (ryx<r.ryx) return -1; else if (ryx>r.ryx) return 1;
  else if (rxz<r.rxz) return -1; else if (rxz>r.rxz) return 1;
  else if (rxy<r.rxy) return -1; else if (rxy>r.rxy) return 1;
  else if (rxx<r.rxx) return -1; else if (rxx>r.rxx) return 1;
  else return 0;
}


const HepRotation HepRotation::IDENTITY;

}  // namespace CLHEP


