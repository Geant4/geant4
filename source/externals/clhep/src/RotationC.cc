// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of methods of the HepRotation class which
// were introduced when ZOOM PhysicsVectors was merged in, which involve
// correcting user-supplied data which is supposed to form a Rotation, or
// rectifying a rotation matrix which may have drifted due to roundoff.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/Rotation.h"

#include <cmath>

namespace CLHEP  {

// --------- Helper methods (private) for setting from 3 columns:

bool HepRotation::setCols 
    ( const Hep3Vector & u1, const Hep3Vector & u2, const Hep3Vector & u3,
      double u1u2,
      Hep3Vector & v1, Hep3Vector & v2, Hep3Vector & v3 ) const {

  if ( (1-std::fabs(u1u2)) <= Hep4RotationInterface::tolerance ) {
    std::cerr << "HepRotation::setCols() - "
      << "All three cols supplied for a Rotation are parallel --"
      << "\n    an arbitrary rotation will be returned" << std::endl;
    setArbitrarily (u1, v1, v2, v3);
    return true;
  }

  v1 = u1;
  v2  = Hep3Vector(u2 - u1u2 * u1).unit();
  v3 = v1.cross(v2);
  if ( v3.dot(u3) >= 0 ) {
    return true;
  } else {
    return false;	// looks more like a reflection in this case!
  }

} // HepRotation::setCols 

void HepRotation::setArbitrarily (const Hep3Vector & ccolX, 
   Hep3Vector & v1, Hep3Vector & v2, Hep3Vector & v3) const {

  // We have all three col's parallel.  Warnings already been given;
  // this just supplies a result which is a valid rotation.

  v1 = ccolX.unit();
  v2 = v1.cross(Hep3Vector(0,0,1));
  if (v2.mag2() != 0) {
    v2 = v2.unit();
  } else {
    v2 = Hep3Vector(1,0,0);
  }
  v3 = v1.cross(v2);

  return;

} // HepRotation::setArbitrarily 



// ----------  Constructors and Assignment:

// 3 orthogonal columns or rows

HepRotation & HepRotation::set( const Hep3Vector & ccolX,
                            	const Hep3Vector & ccolY,
                          	const Hep3Vector & ccolZ ) {
  Hep3Vector ucolX = ccolX.unit();
  Hep3Vector ucolY = ccolY.unit();
  Hep3Vector ucolZ = ccolZ.unit();

  double u1u2 = ucolX.dot(ucolY);
  double f12  = std::fabs(u1u2);
  if ( f12 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepRotation::set() - "
      << "col's X and Y supplied for Rotation are not close to orthogonal"
      << std::endl;
  }
  double u1u3 = ucolX.dot(ucolZ);
  double f13  = std::fabs(u1u3);
  if ( f13 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepRotation::set() - "
      << "col's X and Z supplied for Rotation are not close to orthogonal"
      << std::endl;
  }
  double u2u3 = ucolY.dot(ucolZ);
  double f23  = std::fabs(u2u3);
  if ( f23 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepRotation::set() - "
      << "col's Y and Z supplied for Rotation are not close to orthogonal"
      << std::endl;
  }

  Hep3Vector v1, v2, v3;
  bool isRotation;
  if ( (f12 <= f13) && (f12 <= f23) ) {
    isRotation = setCols ( ucolX, ucolY, ucolZ, u1u2, v1, v2, v3 );
    if ( !isRotation ) {
      std::cerr << "HepRotation::set() - "
        << "col's X Y and Z supplied form closer to a reflection than a Rotation "
        << "\n     col Z is set to col X cross col Y" << std::endl;
    }
  } else if ( f13 <= f23 ) {
    isRotation = setCols ( ucolZ, ucolX, ucolY, u1u3, v3, v1, v2 );
    if ( !isRotation ) {
      std::cerr << "HepRotation::set() - "
        << "col's X Y and Z supplied form closer to a reflection than a Rotation "
        << "\n     col Y is set to col Z cross col X" << std::endl;
    }
  } else {
    isRotation = setCols ( ucolY, ucolZ, ucolX, u2u3, v2, v3, v1 );
    if ( !isRotation ) {
      std::cerr << "HepRotation::set() - "
        << "col's X Y and Z supplied form closer to a reflection than a Rotation "
        << "\n     col X is set to col Y cross col Z" << std::endl;
    }
  }

  rxx = v1.x();  ryx = v1.y(); rzx = v1.z();
  rxy = v2.x();  ryy = v2.y(); rzy = v2.z();
  rxz = v3.x();  ryz = v3.y(); rzz = v3.z();

  return *this;

}  // HepRotation::set(colX, colY, colZ)

HepRotation::HepRotation ( const Hep3Vector & ccolX,
              		   const Hep3Vector & ccolY,
		           const Hep3Vector & ccolZ ) 
{
  set (ccolX, ccolY, ccolZ);
}

HepRotation & HepRotation::setRows( const Hep3Vector & rrowX,
                           	    const Hep3Vector & rrowY,
                              	    const Hep3Vector & rrowZ ) {
  set (rrowX, rrowY, rrowZ);
  invert();
  return *this;
}

// ------- Rectify a near-rotation

void HepRotation::rectify() {
  // Assuming the representation of this is close to a true Rotation,
  // but may have drifted due to round-off error from many operations,
  // this forms an "exact" orthonormal matrix for the rotation again.

  // The first step is to average with the transposed inverse.  This
  // will correct for small errors such as those occuring when decomposing
  // a LorentzTransformation.  Then we take the bull by the horns and
  // formally extract the axis and delta (assuming the Rotation were true)
  // and re-setting the rotation according to those.

  double det =  rxx * ryy * rzz +
                   rxy * ryz * rzx +
                   rxz * ryx * rzy -
                   rxx * ryz * rzy -
                   rxy * ryx * rzz -
                   rxz * ryy * rzx   ;
  if (det <= 0) {
    std::cerr << "HepRotation::rectify() - "
        << "Attempt to rectify a Rotation with determinant <= 0" << std::endl;
    return;
  }
  double di = 1.0 / det;

  // xx, xy, ... are components of inverse matrix:
  double xx1 = (ryy * rzz - ryz * rzy) * di;
  double xy1 = (rzy * rxz - rzz * rxy) * di;
  double xz1 = (rxy * ryz - rxz * ryy) * di;
  double yx1 = (ryz * rzx - ryx * rzz) * di;
  double yy1 = (rzz * rxx - rzx * rxz) * di;
  double yz1 = (rxz * ryx - rxx * ryz) * di;
  double zx1 = (ryx * rzy - ryy * rzx) * di;
  double zy1 = (rzx * rxy - rzy * rxx) * di;
  double zz1 = (rxx * ryy - rxy * ryx) * di;

  // Now average with the TRANSPOSE of that:
  rxx = .5*(rxx + xx1);
  rxy = .5*(rxy + yx1);
  rxz = .5*(rxz + zx1);
  ryx = .5*(ryx + xy1);
  ryy = .5*(ryy + yy1);
  ryz = .5*(ryz + zy1);
  rzx = .5*(rzx + xz1);
  rzy = .5*(rzy + yz1);
  rzz = .5*(rzz + zz1);

  // Now force feed this improved rotation
  double del = delta();
  Hep3Vector u = axis();
  u = u.unit(); // Because if the rotation is inexact, then the
                // axis() returned will not have length 1!
  set(u, del);

} // rectify()

}  // namespace CLHEP

