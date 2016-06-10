// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the HepBoost class.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/Boost.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/LorentzRotation.h"

namespace CLHEP  {

// ----------  Constructors and Assignment:

HepBoost & HepBoost::set (double bx, double by, double bz) {
  double bp2 = bx*bx + by*by + bz*bz;
//  if (bp2 >= 1) {
//    std::cerr << "HepBoost::set() - "
//      << "Boost Vector supplied to set HepBoost represents speed >= c." << std::endl;
//  }    
  double ggamma = 1.0 / std::sqrt(1.0 - bp2);
  double bgamma = ggamma * ggamma / (1.0 + ggamma);
  rep_.xx_ = 1.0 + bgamma * bx * bx;
  rep_.yy_ = 1.0 + bgamma * by * by;
  rep_.zz_ = 1.0 + bgamma * bz * bz;
  rep_.xy_ = bgamma * bx * by;
  rep_.xz_ = bgamma * bx * bz;
  rep_.yz_ = bgamma * by * bz;
  rep_.xt_ = ggamma * bx;
  rep_.yt_ = ggamma * by;
  rep_.zt_ = ggamma * bz;
  rep_.tt_ = ggamma;
  return *this;
}

HepBoost & HepBoost::set (const HepRep4x4Symmetric & m1) {
  rep_ = m1;
  return *this;
}

HepBoost & HepBoost::set (Hep3Vector ddirection, double bbeta) {
  double length = ddirection.mag();
  if (length <= 0) {				// Nan-proofing
    std::cerr << "HepBoost::set() - "
      << "Direction supplied to set HepBoost is zero." << std::endl;
    set (0,0,0);
    return *this;
  }    
  set(bbeta*ddirection.x()/length,
      bbeta*ddirection.y()/length,
      bbeta*ddirection.z()/length);
  return *this;
}

HepBoost & HepBoost::set (const Hep3Vector & bboost) {
  return set (bboost.x(), bboost.y(), bboost.z());
}

// ----------  Accessors:

// ----------  Decomposition:

void HepBoost::decompose (HepRotation & rotation, HepBoost & boost) const {
  HepAxisAngle vdelta = HepAxisAngle();
  rotation = HepRotation(vdelta);
  Hep3Vector bbeta = boostVector();
  boost = HepBoost(bbeta);
}

void HepBoost::decompose (HepAxisAngle & rotation, Hep3Vector & boost) const {
  rotation = HepAxisAngle();
  boost = boostVector();
}

void HepBoost::decompose (HepBoost & boost, HepRotation & rotation) const {
  HepAxisAngle vdelta = HepAxisAngle();
  rotation = HepRotation(vdelta);
  Hep3Vector bbeta = boostVector();
  boost = HepBoost(bbeta);
}

void HepBoost::decompose (Hep3Vector & boost, HepAxisAngle & rotation) const {
  rotation = HepAxisAngle();
  boost = boostVector();
}

// ----------  Comparisons:

double HepBoost::distance2( const HepRotation & r ) const {
  double db2 = norm2();
  double dr2  = r.norm2();
  return (db2 + dr2);
}

double HepBoost::distance2( const HepLorentzRotation & lt ) const {
  HepBoost b1;
  HepRotation r1;
  lt.decompose(b1,r1);
  double db2 = distance2(b1);
  double dr2  = r1.norm2();
  return (db2 + dr2);
}

double HepBoost::howNear ( const HepRotation & r  ) const {
  return std::sqrt(distance2(r));
}

double HepBoost::howNear ( const HepLorentzRotation & lt  ) const {
  return std::sqrt(distance2(lt));
}

bool HepBoost::isNear (const HepRotation & r, double epsilon) const {
  double db2 = norm2();
  if (db2 > epsilon*epsilon) return false;
  double dr2  = r.norm2();
  return (db2+dr2 <= epsilon*epsilon);
}

bool HepBoost::isNear (const HepLorentzRotation & lt, 
			           double epsilon) const {
  HepBoost b1;
  HepRotation r1;
  double db2 = distance2(b1);
  lt.decompose(b1,r1);
  if (db2 > epsilon*epsilon) return false;
  double dr2  = r1.norm2();
  return (db2 + dr2);
}

// ----------  Properties:

double HepBoost::norm2() const {
  double bgx = rep_.xt_;
  double bgy = rep_.yt_;
  double bgz = rep_.zt_;
  return bgx*bgx+bgy*bgy+bgz*bgz;
}

void HepBoost::rectify() {
  // Assuming the representation of this is close to a true pure boost,
  // but may have drifted due to round-off error from many operations,
  // this forms an "exact" pure boost matrix for the LT again.

  // The natural way to do this is to use the t column as a boost and set 
  // based on that boost vector.
  
  // There is perhaps danger that this boost vector will appear equal to or 
  // greater than a unit vector; the best we can do for such a case is use
  // a boost in that direction but rescaled to just less than one.

  // There is in principle no way that gamma could have become negative,
  // but if that happens, we ZMthrow and (if continuing) just rescale, which
  // will change the sign of the last column when computing the boost.

  double gam = tt();
  if (gam <= 0) {				    // 4/12/01 mf 
    std::cerr << "HepBoost::rectify() - "
      << "Attempt to rectify a boost with non-positive gamma." << std::endl;
    if (gam==0) return;				    // NaN-proofing
  }    
  Hep3Vector boost (xt(), yt(), zt());
  boost /= tt();
  if ( boost.mag2() >= 1 ) {			    // NaN-proofing:
    boost /= ( boost.mag() * ( 1.0 + 1.0e-16 ) );   // used to just check > 1
  }
  set ( boost );
}

// ---------- Application is all in .icc 

// ---------- Operations in the group of 4-Rotations

HepLorentzRotation
HepBoost::matrixMultiplication(const HepRep4x4 & m1) const {
  HepRep4x4Symmetric r = rep4x4Symmetric();
  return HepLorentzRotation( HepRep4x4 (
    r.xx_*m1.xx_ + r.xy_*m1.yx_ + r.xz_*m1.zx_ + r.xt_*m1.tx_,
    r.xx_*m1.xy_ + r.xy_*m1.yy_ + r.xz_*m1.zy_ + r.xt_*m1.ty_,
    r.xx_*m1.xz_ + r.xy_*m1.yz_ + r.xz_*m1.zz_ + r.xt_*m1.tz_,
    r.xx_*m1.xt_ + r.xy_*m1.yt_ + r.xz_*m1.zt_ + r.xt_*m1.tt_,

    r.xy_*m1.xx_ + r.yy_*m1.yx_ + r.yz_*m1.zx_ + r.yt_*m1.tx_,
    r.xy_*m1.xy_ + r.yy_*m1.yy_ + r.yz_*m1.zy_ + r.yt_*m1.ty_,
    r.xy_*m1.xz_ + r.yy_*m1.yz_ + r.yz_*m1.zz_ + r.yt_*m1.tz_,
    r.xy_*m1.xt_ + r.yy_*m1.yt_ + r.yz_*m1.zt_ + r.yt_*m1.tt_,

    r.xz_*m1.xx_ + r.yz_*m1.yx_ + r.zz_*m1.zx_ + r.zt_*m1.tx_,
    r.xz_*m1.xy_ + r.yz_*m1.yy_ + r.zz_*m1.zy_ + r.zt_*m1.ty_,
    r.xz_*m1.xz_ + r.yz_*m1.yz_ + r.zz_*m1.zz_ + r.zt_*m1.tz_,
    r.xz_*m1.xt_ + r.yz_*m1.yt_ + r.zz_*m1.zt_ + r.zt_*m1.tt_,

    r.xt_*m1.xx_ + r.yt_*m1.yx_ + r.zt_*m1.zx_ + r.tt_*m1.tx_,
    r.xt_*m1.xy_ + r.yt_*m1.yy_ + r.zt_*m1.zy_ + r.tt_*m1.ty_,
    r.xt_*m1.xz_ + r.yt_*m1.yz_ + r.zt_*m1.zz_ + r.tt_*m1.tz_,
    r.xt_*m1.xt_ + r.yt_*m1.yt_ + r.zt_*m1.zt_ + r.tt_*m1.tt_) );
}

HepLorentzRotation
HepBoost::matrixMultiplication(const HepRep4x4Symmetric & m1) const {
  HepRep4x4Symmetric r = rep4x4Symmetric();
  return HepLorentzRotation( HepRep4x4 (
    r.xx_*m1.xx_ + r.xy_*m1.xy_ + r.xz_*m1.xz_ + r.xt_*m1.xt_,
    r.xx_*m1.xy_ + r.xy_*m1.yy_ + r.xz_*m1.yz_ + r.xt_*m1.yt_,
    r.xx_*m1.xz_ + r.xy_*m1.yz_ + r.xz_*m1.zz_ + r.xt_*m1.zt_,
    r.xx_*m1.xt_ + r.xy_*m1.yt_ + r.xz_*m1.zt_ + r.xt_*m1.tt_,

    r.xy_*m1.xx_ + r.yy_*m1.xy_ + r.yz_*m1.xz_ + r.yt_*m1.xt_,
    r.xy_*m1.xy_ + r.yy_*m1.yy_ + r.yz_*m1.yz_ + r.yt_*m1.yt_,
    r.xy_*m1.xz_ + r.yy_*m1.yz_ + r.yz_*m1.zz_ + r.yt_*m1.zt_,
    r.xy_*m1.xt_ + r.yy_*m1.yt_ + r.yz_*m1.zt_ + r.yt_*m1.tt_,

    r.xz_*m1.xx_ + r.yz_*m1.xy_ + r.zz_*m1.xz_ + r.zt_*m1.xt_,
    r.xz_*m1.xy_ + r.yz_*m1.yy_ + r.zz_*m1.yz_ + r.zt_*m1.yt_,
    r.xz_*m1.xz_ + r.yz_*m1.yz_ + r.zz_*m1.zz_ + r.zt_*m1.zt_,
    r.xz_*m1.xt_ + r.yz_*m1.yt_ + r.zz_*m1.zt_ + r.zt_*m1.tt_,

    r.xt_*m1.xx_ + r.yt_*m1.xy_ + r.zt_*m1.xz_ + r.tt_*m1.xt_,
    r.xt_*m1.xy_ + r.yt_*m1.yy_ + r.zt_*m1.yz_ + r.tt_*m1.yt_,
    r.xt_*m1.xz_ + r.yt_*m1.yz_ + r.zt_*m1.zz_ + r.tt_*m1.zt_,
    r.xt_*m1.xt_ + r.yt_*m1.yt_ + r.zt_*m1.zt_ + r.tt_*m1.tt_) );
}

HepLorentzRotation
HepBoost::operator* (const HepLorentzRotation & lt) const {
  return matrixMultiplication(lt.rep4x4());
}

HepLorentzRotation
HepBoost::operator* (const HepBoost & b) const {
  return matrixMultiplication(b.rep_);
}

HepLorentzRotation
HepBoost::operator* (const HepRotation & r) const {
  return matrixMultiplication(r.rep4x4());
}

// ---------- I/O:

std::ostream & HepBoost::print( std::ostream & os ) const {
  if ( rep_.tt_ <= 1 ) {
    os << "Lorentz Boost( IDENTITY )";
  } else {
    double norm = boostVector().mag();
    os << "\nLorentz Boost " << boostVector()/norm <<
          "\n{beta = " << beta() << " gamma = " << gamma() << "}\n";
  }
  return os;
}

}  // namespace CLHEP
