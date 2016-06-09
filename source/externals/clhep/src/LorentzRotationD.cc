// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of those parts of the HepLorentzRotation class
// which involve decomposition into Boost*Rotation.

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/LorentzRotation.h"

namespace CLHEP  {

// ----------  Decomposition:

void HepLorentzRotation::decompose 
	(HepBoost & bboost, HepRotation & rotation) const {

  // The boost will be the pure boost based on column 4 of the transformation
  // matrix.  Since the constructor takes the beta vector, and not beta*gamma,
  // we first divide through by gamma = the tt element.  This of course can
  // never be zero since the last row has t**2 - v**2 = +1.

  Hep3Vector betaVec ( xt(), yt(), zt() );
  betaVec *= 1.0 / tt();
  bboost.set( betaVec );

  // The rotation will be inverse of B times T.

  HepBoost B( -betaVec );
  HepLorentzRotation R( B * *this );

  HepRep3x3 m1  ( R.xx(), R.xy(), R.xz(),
                  R.yx(), R.yy(), R.yz(),
                  R.zx(), R.zy(), R.zz() );
  rotation.set( m1 );
  rotation.rectify();
  
  return;

}

void HepLorentzRotation::decompose 
	(Hep3Vector & bboost, HepAxisAngle & rotation) const {
  HepRotation r;
  HepBoost b;
  decompose(b,r);
  bboost = b.boostVector();
  rotation = r.axisAngle();
  return;
}

void HepLorentzRotation::decompose 
	(HepRotation & rotation, HepBoost & bboost) const {

  // In this case the pure boost is based on row 4 of the matrix.  

  Hep3Vector betaVec( tx(), ty(), tz() );
  betaVec *= 1.0 / tt();
  bboost.set( betaVec );

  // The rotation will be T times the inverse of B.

  HepBoost B( -betaVec );
  HepLorentzRotation R( *this * B );

  HepRep3x3 m1 ( R.xx(), R.xy(), R.xz(),
                 R.yx(), R.yy(), R.yz(),
                 R.zx(), R.zy(), R.zz() );
  rotation.set( m1 );
  rotation.rectify();
  return;

}

void HepLorentzRotation::decompose 
	(HepAxisAngle & rotation, Hep3Vector & bboost) const {
  HepRotation r;
  HepBoost b;
  decompose(r,b);
  rotation = r.axisAngle();
  bboost = b.boostVector();
  return;
}

double HepLorentzRotation::distance2( const HepBoost & b ) const {
  HepBoost    b1;
  HepRotation r1; 
  decompose( b1, r1 );
  double db2 = b1.distance2( b );
  double dr2 = r1.norm2(); 
  return ( db2 + dr2 );
}

double HepLorentzRotation::distance2( const HepRotation & r ) const {
  HepBoost    b1;
  HepRotation r1; 
  decompose( b1, r1 );
  double db2 = b1.norm2( );
  double dr2 = r1.distance2( r ); 
  return ( db2 + dr2 );
}

double HepLorentzRotation::distance2( 
				   const HepLorentzRotation & lt  ) const {
  HepBoost    b1;
  HepRotation r1; 
  decompose( b1, r1 );
  HepBoost    b2;
  HepRotation r2; 
  lt.decompose (b2, r2);
  double db2 = b1.distance2( b2 );
  double dr2 = r1.distance2( r2 ); 
  return ( db2 + dr2 );
}

double HepLorentzRotation::howNear( const HepBoost & b ) const {
  return std::sqrt( distance2( b ) );
}
double HepLorentzRotation::howNear( const HepRotation & r ) const {
  return std::sqrt( distance2( r ) );
}
double HepLorentzRotation::howNear( const HepLorentzRotation & lt )const {
  return std::sqrt( distance2( lt ) );
}

bool HepLorentzRotation::isNear(
		const HepBoost & b, double epsilon ) const {
  HepBoost    b1;
  HepRotation r1; 
  decompose( b1, r1 );
  double db2 = b1.distance2(b);
  if ( db2 > epsilon*epsilon ) {
     return false;       // Saves the time-consuming Rotation::norm2
  }
  double dr2 = r1.norm2();
  return ( (db2 + dr2) <= epsilon*epsilon );
}

bool HepLorentzRotation::isNear(
		const HepRotation & r, double epsilon ) const {
  HepBoost    b1;
  HepRotation r1; 
  decompose( b1, r1 );
  double db2 = b1.norm2();
  if ( db2 > epsilon*epsilon ) {
     return false;       // Saves the time-consuming Rotation::distance2
  }
  double dr2 = r1.distance2(r);
  return ( (db2 + dr2) <= epsilon*epsilon );
}

bool HepLorentzRotation::isNear(
		const HepLorentzRotation & lt, double epsilon ) const {
  HepBoost    b1;
  HepRotation r1; 
  decompose( b1, r1 );
  HepBoost    b2;
  HepRotation r2; 
  lt.decompose (b2, r2);
  double db2 = b1.distance2(b2);
  if ( db2 > epsilon*epsilon ) {
     return false;       // Saves the time-consuming Rotation::distance2
  }
  double dr2 = r1.distance2(r2);
  return ( (db2 + dr2) <= epsilon*epsilon );
}

double HepLorentzRotation::norm2() const {
  HepBoost    b;
  HepRotation r;
  decompose( b, r );
  return b.norm2() + r.norm2();
}

void HepLorentzRotation::rectify() {
  
  // Assuming the representation of this is close to a true LT,
  // but may have drifted due to round-off error from many operations,
  // this forms an "exact" orthosymplectic matrix for the LT again.
 
  // There are several ways to do this, all equivalent to lowest order in
  // the corrected error.  We choose to form an LT based on the inverse boost
  // extracted from row 4, and left-multiply by LT to form what would be
  // a rotation if the LT were kosher.  We drop the possibly non-zero t
  // components of that, rectify that rotation and multiply back by the boost.
                    
  Hep3Vector beta (tx(), ty(), tz());
  double gam = tt();			// NaN-proofing
  if ( gam <= 0 ) {
    std::cerr << "HepLorentzRotation::rectify() - "
	<< "rectify() on a transformation with tt() <= 0 - will not help!"
        << std::endl;
    gam = 1;
  }
  beta *= 1.0/gam;
  HepLorentzRotation R = (*this) * HepBoost(-beta);

  HepRep3x3  m1 ( R.xx(), R.xy(), R.xz(),
                  R.yx(), R.yy(), R.yz(),
                  R.zx(), R.zy(), R.zz() );

  HepRotation Rgood (m1);
  Rgood.rectify();

  set ( Rgood, HepBoost(beta) );
}

}  // namespace CLHEP
