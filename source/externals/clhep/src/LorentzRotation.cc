// -*- C++ -*-
// $Id:$
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation basic parts of the HepLorentzRotation class.
//
// Some ZOOM methods involving construction from columns and decomposition 
// into boost*rotation are split off into LorentzRotationC and LorentzRotationD

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/LorentzRotation.h"

#include <iostream>
#include <iomanip>

namespace CLHEP  {

// ----------  Constructors and Assignment:


HepLorentzRotation & HepLorentzRotation::set
				(double bx, double by, double bz) {
  double bp2 = bx*bx + by*by + bz*bz;
//  if (bp2 >= 1) {
//    std::cerr << "HepLorentzRotation::set() - "
//      << "Boost Vector supplied to set HepLorentzRotation represents speed >= c."
//      << std::endl;
//  }    
  double gamma = 1.0 / std::sqrt(1.0 - bp2);
  double bgamma = gamma * gamma / (1.0 + gamma);
  mxx = 1.0 + bgamma * bx * bx;
  myy = 1.0 + bgamma * by * by;
  mzz = 1.0 + bgamma * bz * bz;
  mxy = myx = bgamma * bx * by;
  mxz = mzx = bgamma * bx * bz;
  myz = mzy = bgamma * by * bz;
  mxt = mtx = gamma * bx;
  myt = mty = gamma * by;
  mzt = mtz = gamma * bz;
  mtt = gamma;
  return *this;
}

HepLorentzRotation & HepLorentzRotation::set 
		(const HepBoost & B, const HepRotation & R) {
  set (B.rep4x4());
  *this = matrixMultiplication ( R.rep4x4() );
  return *this;
}

HepLorentzRotation & HepLorentzRotation::set 
		(const HepRotation & R, const HepBoost & B) {
  set (R.rep4x4());
  *this = matrixMultiplication ( B.rep4x4() );
  return *this;
}

// ----------  Accessors:

// ------------  Subscripting:

double HepLorentzRotation::operator () (int i, int j) const {
  if (i == 0) {
    if (j == 0) { return xx(); }
    if (j == 1) { return xy(); }
    if (j == 2) { return xz(); } 
    if (j == 3) { return xt(); } 
  } else if (i == 1) {
    if (j == 0) { return yx(); }
    if (j == 1) { return yy(); }
    if (j == 2) { return yz(); } 
    if (j == 3) { return yt(); } 
  } else if (i == 2) {
    if (j == 0) { return zx(); }
    if (j == 1) { return zy(); }
    if (j == 2) { return zz(); } 
    if (j == 3) { return zt(); } 
  } else if (i == 3) {
    if (j == 0) { return tx(); }
    if (j == 1) { return ty(); }
    if (j == 2) { return tz(); } 
    if (j == 3) { return tt(); } 
  } 
  std::cerr << "HepLorentzRotation subscripting: bad indeces "
	    << "(" << i << "," << j << ")\n";
  return 0.0;
} 

// ---------- Application:


// ---------- Comparison:

int HepLorentzRotation::compare( const HepLorentzRotation & m1  ) const {
       if (mtt<m1.mtt) return -1; else if (mtt>m1.mtt) return 1;
  else if (mtz<m1.mtz) return -1; else if (mtz>m1.mtz) return 1;
  else if (mty<m1.mty) return -1; else if (mty>m1.mty) return 1;
  else if (mtx<m1.mtx) return -1; else if (mtx>m1.mtx) return 1;

  else if (mzt<m1.mzt) return -1; else if (mzt>m1.mzt) return 1;
  else if (mzz<m1.mzz) return -1; else if (mzz>m1.mzz) return 1;
  else if (mzy<m1.mzy) return -1; else if (mzy>m1.mzy) return 1;
  else if (mzx<m1.mzx) return -1; else if (mzx>m1.mzx) return 1;

  else if (myt<m1.myt) return -1; else if (myt>m1.myt) return 1;
  else if (myz<m1.myz) return -1; else if (myz>m1.myz) return 1;
  else if (myy<m1.myy) return -1; else if (myy>m1.myy) return 1;
  else if (myx<m1.myx) return -1; else if (myx>m1.myx) return 1;

  else if (mxt<m1.mxt) return -1; else if (mxt>m1.mxt) return 1;
  else if (mxz<m1.mxz) return -1; else if (mxz>m1.mxz) return 1;
  else if (mxy<m1.mxy) return -1; else if (mxy>m1.mxy) return 1;
  else if (mxx<m1.mxx) return -1; else if (mxx>m1.mxx) return 1;

  else return 0;
}


// ---------- Operations in the group of 4-Rotations

HepLorentzRotation
HepLorentzRotation::matrixMultiplication(const HepRep4x4 & m1) const {
  return HepLorentzRotation(
    mxx*m1.xx_ + mxy*m1.yx_ + mxz*m1.zx_ + mxt*m1.tx_,
    mxx*m1.xy_ + mxy*m1.yy_ + mxz*m1.zy_ + mxt*m1.ty_,
    mxx*m1.xz_ + mxy*m1.yz_ + mxz*m1.zz_ + mxt*m1.tz_,
    mxx*m1.xt_ + mxy*m1.yt_ + mxz*m1.zt_ + mxt*m1.tt_,

    myx*m1.xx_ + myy*m1.yx_ + myz*m1.zx_ + myt*m1.tx_,
    myx*m1.xy_ + myy*m1.yy_ + myz*m1.zy_ + myt*m1.ty_,
    myx*m1.xz_ + myy*m1.yz_ + myz*m1.zz_ + myt*m1.tz_,
    myx*m1.xt_ + myy*m1.yt_ + myz*m1.zt_ + myt*m1.tt_,

    mzx*m1.xx_ + mzy*m1.yx_ + mzz*m1.zx_ + mzt*m1.tx_,
    mzx*m1.xy_ + mzy*m1.yy_ + mzz*m1.zy_ + mzt*m1.ty_,
    mzx*m1.xz_ + mzy*m1.yz_ + mzz*m1.zz_ + mzt*m1.tz_,
    mzx*m1.xt_ + mzy*m1.yt_ + mzz*m1.zt_ + mzt*m1.tt_,

    mtx*m1.xx_ + mty*m1.yx_ + mtz*m1.zx_ + mtt*m1.tx_,
    mtx*m1.xy_ + mty*m1.yy_ + mtz*m1.zy_ + mtt*m1.ty_,
    mtx*m1.xz_ + mty*m1.yz_ + mtz*m1.zz_ + mtt*m1.tz_,
    mtx*m1.xt_ + mty*m1.yt_ + mtz*m1.zt_ + mtt*m1.tt_ );
}

HepLorentzRotation & HepLorentzRotation::rotateX(double delta) {
  double c1 = std::cos (delta);
  double s1 = std::sin (delta);
  HepLorentzVector rowy = row2();
  HepLorentzVector rowz = row3();
  HepLorentzVector r2 = c1 * rowy - s1 * rowz;
  HepLorentzVector r3 = s1 * rowy + c1 * rowz;
  myx = r2.x();   myy = r2.y();   myz = r2.z();   myt = r2.t();	
  mzx = r3.x();   mzy = r3.y();   mzz = r3.z();   mzt = r3.t();	
  return *this;
}

HepLorentzRotation & HepLorentzRotation::rotateY(double delta) {
  double c1 = std::cos (delta);
  double s1 = std::sin (delta);
  HepLorentzVector rowx = row1();
  HepLorentzVector rowz = row3();
  HepLorentzVector r1 =  c1 * rowx + s1 * rowz;
  HepLorentzVector r3 = -s1 * rowx + c1 * rowz;
  mxx = r1.x();   mxy = r1.y();   mxz = r1.z();   mxt = r1.t();	
  mzx = r3.x();   mzy = r3.y();   mzz = r3.z();   mzt = r3.t();	
  return *this;
}

HepLorentzRotation & HepLorentzRotation::rotateZ(double delta) {
  double c1 = std::cos (delta);
  double s1 = std::sin (delta);
  HepLorentzVector rowx = row1();
  HepLorentzVector rowy = row2();
  HepLorentzVector r1 = c1 * rowx - s1 * rowy;
  HepLorentzVector r2 = s1 * rowx + c1 * rowy;
  mxx = r1.x();   mxy = r1.y();   mxz = r1.z();   mxt = r1.t();
  myx = r2.x();   myy = r2.y();   myz = r2.z();   myt = r2.t();
  return *this;
}

HepLorentzRotation & HepLorentzRotation::boostX(double beta) {
  double b2 = beta*beta;
//  if (b2 >= 1) {
//    std::cerr << "HepLorentzRotation::boostX() - "
//      << "Beta supplied to HepLorentzRotation::boostX represents speed >= c."
//      << std::endl;
//  }    
  double g1  = 1.0/std::sqrt(1.0-b2);
  double bg = beta*g1;
  HepLorentzVector rowx = row1();
  HepLorentzVector rowt = row4();
  HepLorentzVector r1 =  g1 * rowx + bg * rowt;
  HepLorentzVector r4 = bg * rowx +  g1 * rowt;
  mxx = r1.x();   mxy = r1.y();   mxz = r1.z();   mxt = r1.t();	
  mtx = r4.x();   mty = r4.y();   mtz = r4.z();   mtt = r4.t();	
  return *this;
}

HepLorentzRotation & HepLorentzRotation::boostY(double beta) {
  double b2 = beta*beta;
//  if (b2 >= 1) {
//    std::cerr << "HepLorentzRotation::boostY() - "
//      << "Beta supplied to HepLorentzRotation::boostY represents speed >= c."
//      << std::endl;
//  }    
  double g1  = 1.0/std::sqrt(1.0-b2);
  double bg = beta*g1;
  HepLorentzVector rowy = row2();
  HepLorentzVector rowt = row4();
  HepLorentzVector r2 =  g1 * rowy + bg * rowt;
  HepLorentzVector r4 = bg * rowy +  g1 * rowt;
  myx = r2.x();   myy = r2.y();   myz = r2.z();   myt = r2.t();	
  mtx = r4.x();   mty = r4.y();   mtz = r4.z();   mtt = r4.t();	
  return *this;
}

HepLorentzRotation & HepLorentzRotation::boostZ(double beta) {
  double b2 = beta*beta;
//  if (b2 >= 1) {
//    std::cerr << "HepLorentzRotation::boostZ() - "
//      << "Beta supplied to HepLorentzRotation::boostZ represents speed >= c."
//      << std::endl;
//  }    
  double g1  = 1.0/std::sqrt(1.0-b2);
  double bg = beta*g1;
  HepLorentzVector rowz = row3();
  HepLorentzVector rowt = row4();
  HepLorentzVector r3 =  g1 * rowz + bg * rowt;
  HepLorentzVector r4 = bg * rowz +  g1 * rowt;
  mtx = r4.x();   mty = r4.y();   mtz = r4.z();   mtt = r4.t();	
  mzx = r3.x();   mzy = r3.y();   mzz = r3.z();   mzt = r3.t();	
  return *this;
}

std::ostream & HepLorentzRotation::print( std::ostream & os ) const {
  os << "\n   [ ( " <<
        std::setw(11) << std::setprecision(6) << xx() << "   " <<
        std::setw(11) << std::setprecision(6) << xy() << "   " <<
        std::setw(11) << std::setprecision(6) << xz() << "   " <<
        std::setw(11) << std::setprecision(6) << xt() << ")\n"
     << "     ( " <<
        std::setw(11) << std::setprecision(6) << yx() << "   " <<
        std::setw(11) << std::setprecision(6) << yy() << "   " <<
        std::setw(11) << std::setprecision(6) << yz() << "   " <<
        std::setw(11) << std::setprecision(6) << yt() << ")\n"
     << "     ( " <<
        std::setw(11) << std::setprecision(6) << zx() << "   " <<
        std::setw(11) << std::setprecision(6) << zy() << "   " <<
        std::setw(11) << std::setprecision(6) << zz() << "   " <<
        std::setw(11) << std::setprecision(6) << zt() << ")\n"
     << "     ( " <<
        std::setw(11) << std::setprecision(6) << tx() << "   " <<
        std::setw(11) << std::setprecision(6) << ty() << "   " <<
        std::setw(11) << std::setprecision(6) << tz() << "   " <<
        std::setw(11) << std::setprecision(6) << tt() << ") ]\n";
  return os;
}

HepLorentzRotation operator* ( const HepRotation & r,
                               const HepLorentzRotation & lt) {
  r.rep4x4();
  lt.rep4x4();
  return HepLorentzRotation( HepRep4x4(
         r.xx()*lt.xx() + r.xy()*lt.yx() + r.xz()*lt.zx() + r.xt()*lt.tx(),
	 r.xx()*lt.xy() + r.xy()*lt.yy() + r.xz()*lt.zy() + r.xt()*lt.ty(),
	 r.xx()*lt.xz() + r.xy()*lt.yz() + r.xz()*lt.zz() + r.xt()*lt.tz(),
	 r.xx()*lt.xt() + r.xy()*lt.yt() + r.xz()*lt.zt() + r.xt()*lt.tt(),

         r.yx()*lt.xx() + r.yy()*lt.yx() + r.yz()*lt.zx() + r.yt()*lt.tx(),
         r.yx()*lt.xy() + r.yy()*lt.yy() + r.yz()*lt.zy() + r.yt()*lt.ty(),
         r.yx()*lt.xz() + r.yy()*lt.yz() + r.yz()*lt.zz() + r.yt()*lt.tz(),
         r.yx()*lt.xt() + r.yy()*lt.yt() + r.yz()*lt.zt() + r.yt()*lt.tt(),

         r.zx()*lt.xx() + r.zy()*lt.yx() + r.zz()*lt.zx() + r.zt()*lt.tx(),
         r.zx()*lt.xy() + r.zy()*lt.yy() + r.zz()*lt.zy() + r.zt()*lt.ty(),
         r.zx()*lt.xz() + r.zy()*lt.yz() + r.zz()*lt.zz() + r.zt()*lt.tz(),
         r.zx()*lt.xt() + r.zy()*lt.yt() + r.zz()*lt.zt() + r.zt()*lt.tt(),

         r.tx()*lt.xx() + r.ty()*lt.yx() + r.tz()*lt.zx() + r.tt()*lt.tx(),
         r.tx()*lt.xy() + r.ty()*lt.yy() + r.tz()*lt.zy() + r.tt()*lt.ty(),
         r.tx()*lt.xz() + r.ty()*lt.yz() + r.tz()*lt.zz() + r.tt()*lt.tz(),
         r.tx()*lt.xt() + r.ty()*lt.yt() + r.tz()*lt.zt() + r.tt()*lt.tt() ) );
}


const HepLorentzRotation HepLorentzRotation::IDENTITY;

}  // namespace CLHEP
