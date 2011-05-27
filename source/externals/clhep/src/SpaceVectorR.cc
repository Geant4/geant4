// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the subset of those methods of the Hep3Vector 
// class which originated from the ZOOM SpaceVector class *and* which involve
// the concepts of rotation.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/AxisAngle.h"
#include "CLHEP/Vector/EulerAngles.h"

namespace CLHEP  {

//-************************
// rotate about axis
//-************************

Hep3Vector & Hep3Vector::rotate (const Hep3Vector & axis,
				   double delta) {
  double r = axis.mag();
  if ( r == 0 ) {
    std::cerr << "Hep3Vector::rotate() - "
      << "Attempt to rotate around a zero vector axis! " << std::endl;
    return *this;
  }
  register double scale=1.0/r;
  register double ux = scale*axis.getX();
  register double uy = scale*axis.getY();
  register double uz = scale*axis.getZ();
  double cd = std::cos(delta);
  double sd = std::sin(delta);
  register double ocd = 1 - cd;
  double rx;
  double ry;
  double rz;

  { register double  ocdux = ocd * ux;
    rx = dx * ( cd + ocdux * ux           ) +
         dy * (      ocdux * uy - sd * uz ) +
         dz * (      ocdux * uz + sd * uy ) ;
  }

  { register double  ocduy = ocd * uy;
    ry = dy * ( cd + ocduy * uy           ) +
         dz * (      ocduy * uz - sd * ux ) +
         dx * (      ocduy * ux + sd * uz ) ;
  }

  { register double  ocduz = ocd * uz;
    rz = dz * ( cd + ocduz * uz           ) +
         dx * (      ocduz * ux - sd * uy ) +
         dy * (      ocduz * uy + sd * ux ) ;
  }

  dx = rx;
  dy = ry;
  dz = rz;

  return *this;
} /* rotate */

//-****************************
// rotate by three euler angles
//-****************************


Hep3Vector & Hep3Vector::rotate (double phi, 
				 double theta, 
				 double psi)  {

  double rx;
  double ry;
  double rz;

  register double sinPhi   = std::sin( phi   ), cosPhi   = std::cos( phi   );
  register double sinTheta = std::sin( theta ), cosTheta = std::cos( theta );
  register double sinPsi   = std::sin( psi   ), cosPsi   = std::cos( psi   );

  rx = 	(cosPsi * cosPhi   - cosTheta * sinPsi * sinPhi)   * dx  +
	(cosPsi * sinPhi   + cosTheta * sinPsi * cosPhi)   * dy  +
  	(sinPsi * sinTheta)				   * dz  ;

  ry = 	(- sinPsi * cosPhi - cosTheta * cosPsi * sinPhi)   * dx  +
	(- sinPsi * sinPhi + cosTheta * cosPsi * cosPhi)   * dy  +
  	(cosPsi * sinTheta)				   * dz  ;

  rz = 	(sinTheta * sinPhi)				   * dx  +
  	(- sinTheta * cosPhi)				   * dy  +
	(cosTheta)					   * dz  ;

  dx = rx;
  dy = ry;
  dz = rz;

  return *this;

} /* rotate */


//-*******************
// rotate(HepAxisAngle)
// rotate(HepEulerAngles)
//-*******************

Hep3Vector & Hep3Vector::rotate (const HepAxisAngle & ax ) {
  return rotate( ax.getAxis(), ax.delta() );
}

Hep3Vector & Hep3Vector::rotate (const HepEulerAngles & ex ) {
  return rotate( ex.phi(), ex.theta(), ex.psi() );
}


//-***********************
// rotationOf(HepAxisAngle)
// rotationOf(HepEulerAngles)
// and coordinate axis rotations
//-***********************

Hep3Vector rotationOf (const Hep3Vector & vec, const HepAxisAngle & ax) {
  Hep3Vector vv(vec);
  return vv.rotate (ax);
}

Hep3Vector rotationOf (const Hep3Vector & vec,
                       const Hep3Vector & axis, double delta) {
  Hep3Vector vv(vec);
  return vv.rotate(axis, delta);
}

Hep3Vector rotationOf (const Hep3Vector & vec, const HepEulerAngles & ex) {
  Hep3Vector vv(vec);
  return vv.rotate (ex);
}

Hep3Vector rotationOf (const Hep3Vector & vec,
                       double phi, double theta, double psi) {
  Hep3Vector vv(vec);
  return vv.rotate(phi, theta, psi);
}

Hep3Vector rotationXOf (const Hep3Vector & vec, double delta) {
  Hep3Vector vv(vec);
  return vv.rotateX (delta);
}

Hep3Vector rotationYOf (const Hep3Vector & vec, double delta) {
  Hep3Vector vv(vec);
  return vv.rotateY (delta);
}

Hep3Vector rotationZOf (const Hep3Vector & vec, double delta) {
  Hep3Vector vv(vec);
  return vv.rotateZ (delta);
}

}  // namespace CLHEP
