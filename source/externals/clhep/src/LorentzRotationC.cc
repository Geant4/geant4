// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of that part of the HepLorentzRotation class
// which is concerned with setting or constructing the transformation based 
// on 4 supplied columns or rows.

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/LorentzVector.h"

#include <cmath>

namespace CLHEP  {

// ----------  Constructors and Assignment:

HepLorentzRotation & HepLorentzRotation::set (const HepLorentzVector & ccol1,
                       			      const HepLorentzVector & ccol2,
			                      const HepLorentzVector & ccol3,
			                      const HepLorentzVector & ccol4) {
  // First, test that the four cols do represent something close to a
  // true LT:

  ZMpvMetric_t savedMetric = HepLorentzVector::setMetric (TimePositive);

  if ( ccol4.getT() < 0 ) {
    std::cerr << "HepLorentzRotation::set() - "
      << "column 4 supplied to define transformation has negative T component"
      << std::endl;
    *this = HepLorentzRotation();
    return *this;
  }
/*
  double u1u1 = ccol1.dot(ccol1);
  double f11  = std::fabs(u1u1 + 1.0);
  if ( f11 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepLorentzRotation::set() - "
      << "column 1 supplied for HepLorentzRotation has w*w != -1" << std::endl;
  }
  double u2u2 = ccol2.dot(ccol2);
  double f22  = std::fabs(u2u2 + 1.0);
  if ( f22 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepLorentzRotation::set() - "
      << "column 2 supplied for HepLorentzRotation has w*w != -1" << std::endl;
  }
  double u3u3 = ccol3.dot(ccol3);
  double f33  = std::fabs(u3u3 + 1.0);
  if ( f33 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepLorentzRotation::set() - "
      << "column 3 supplied for HepLorentzRotation has w*w != -1" << std::endl;
  }
  double u4u4 = ccol4.dot(ccol4);
  double f44  = std::fabs(u4u4 - 1.0);
  if ( f44 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepLorentzRotation::set() - "
      << "column 4 supplied for HepLorentzRotation has w*w != +1" << std::endl;
  }

  double u1u2 = ccol1.dot(ccol2);
  double f12  = std::fabs(u1u2);
  if ( f12 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepLorentzRotation::set() - "
      << "columns 1 and 2 supplied for HepLorentzRotation have non-zero dot" << std::endl;
  }
  double u1u3 = ccol1.dot(ccol3);
  double f13  = std::fabs(u1u3);

  if ( f13 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepLorentzRotation::set() - "
      << "columns 1 and 3 supplied for HepLorentzRotation have non-zero dot" << std::endl;
  }
  double u1u4 = ccol1.dot(ccol4);
  double f14  = std::fabs(u1u4);
  if ( f14 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepLorentzRotation::set() - "
      << "columns 1 and 4 supplied for HepLorentzRotation have non-zero dot" << std::endl;
  }
  double u2u3 = ccol2.dot(ccol3);
  double f23  = std::fabs(u2u3);
  if ( f23 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepLorentzRotation::set() - "
      << "columns 2 and 3 supplied for HepLorentzRotation have non-zero dot" << std::endl;
  }
  double u2u4 = ccol2.dot(ccol4);
  double f24  = std::fabs(u2u4);
  if ( f24 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepLorentzRotation::set() - "
      << "columns 2 and 4 supplied for HepLorentzRotation have non-zero dot" << std::endl;
  }
  double u3u4 = ccol3.dot(ccol4);
  double f34  = std::fabs(u3u4);
  if ( f34 > Hep4RotationInterface::tolerance ) {
    std::cerr << "HepLorentzRotation::set() - "
      << "columns 3 and 4 supplied for HepLorentzRotation have non-zero dot" << std::endl;
  }
*/
  // Our strategy will be to order the cols, then do gram-schmidt on them
  // (that is, remove the components of col d that make it non-orthogonal to
  // col c, normalize that, then remove the components of b that make it
  // non-orthogonal to d and to c, normalize that, etc.

  // Because col4, the time col, is most likely to be computed directly, we
  // will start from there and work left-ward.

  HepLorentzVector a, b, c, d;
  bool isLorentzTransformation = true;
  double norm;

  d = ccol4;
  norm = d.dot(d);
  if (norm <= 0.0) {
    isLorentzTransformation = false;
    if (norm == 0.0) {
      d = T_HAT4;       // Moot, but let's keep going...
      norm = 1.0;
    }
  }
  d /= norm;

  c = ccol3 - ccol3.dot(d) * d;
  norm = -c.dot(c);
  if (norm <= 0.0) {
    isLorentzTransformation = false;
    if (norm == 0.0) {
      c = Z_HAT4;       // Moot
      norm = 1.0;
    }
  }
  c /= norm;

  b = ccol2 + ccol2.dot(c) * c - ccol2.dot(d) * d;
  norm = -b.dot(b);
  if (norm <= 0.0) {
    isLorentzTransformation = false;
    if (norm == 0.0) {
      b = Y_HAT4;       // Moot
      norm = 1.0;
    }
  }
  b /= norm;

  a = ccol1 + ccol1.dot(b) * b + ccol1.dot(c) * c - ccol1.dot(d) * d;
  norm = -a.dot(a);
  if (norm <= 0.0) {
    isLorentzTransformation = false;
    if (norm == 0.0) {
      a = X_HAT4;       // Moot
      norm = 1.0;
    }
  }
  a /= norm;

  if ( !isLorentzTransformation ) {
      std::cerr << "HepLorentzRotation::set() - "
        << "cols 1-4 supplied to define transformation form either \n"
        << "       a boosted reflection or a tachyonic transformation -- \n"
        << "       transformation will be set to Identity " << std::endl;

    *this = HepLorentzRotation();
  }

  if ( isLorentzTransformation ) {
    mxx = a.x(); myx = a.y(); mzx = a.z(); mtx = a.t();
    mxy = b.x(); myy = b.y(); mzy = b.z(); mty = b.t();
    mxz = c.x(); myz = c.y(); mzz = c.z(); mtz = c.t();
    mxt = d.x(); myt = d.y(); mzt = d.z(); mtt = d.t();
  }

  HepLorentzVector::setMetric (savedMetric);
  return *this;

} // set ( col1, col2, col3, col4 )

HepLorentzRotation & HepLorentzRotation::setRows
	 (const HepLorentzVector & rrow1,
          const HepLorentzVector & rrow2,
	  const HepLorentzVector & rrow3,
	  const HepLorentzVector & rrow4) {
  // Set based on using those rows as columns:
  set (rrow1, rrow2, rrow3, rrow4);
  // Now transpose in place:
  double q1, q2, q3;
  q1  = mxy;  q2  = mxz;  q3  = mxt;
  mxy = myx;  mxz = mzx;  mxt = mtx;
  myx = q1;   mzx = q2;   mtx = q3;
  q1  = myz;  q2  = myt;  q3  = mzt;
  myz = mzy;  myt = mty;  mzt = mtz;
  mzy = q1;   mty = q2;   mtz = q3;
  return *this;
} // LorentzTransformation::setRows(row1 ... row4)

HepLorentzRotation::HepLorentzRotation ( const HepLorentzVector & ccol1,
		                         const HepLorentzVector & ccol2,
                		         const HepLorentzVector & ccol3,
                       			 const HepLorentzVector & ccol4 ) 
{
  set ( ccol1, ccol2, ccol3, ccol4 );
}

}  // namespace CLHEP

