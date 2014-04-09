// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is part of the implementation of the HepLorentzVector class:
// Those methods which originated from ZOOM and which deal with relativistic
// kinematic properties.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/LorentzVector.h"

#include <cmath>

namespace CLHEP  {

//-******************
// Metric flexibility
//-******************

ZMpvMetric_t HepLorentzVector::setMetric( ZMpvMetric_t a1 ) {
  ZMpvMetric_t oldMetric = (metric > 0) ? TimePositive : TimeNegative;
  if ( a1 == TimeNegative ) {
    metric = -1.0;
  } else {
    metric =  1.0;
  }
  return oldMetric;
}

ZMpvMetric_t HepLorentzVector::getMetric() {
  return ( (metric > 0) ? TimePositive : TimeNegative );
}

//-********
// plus
// minus
//-********

double HepLorentzVector::plus (const Hep3Vector & ref) const {
  double r = ref.mag();
  if (r == 0) {
    std::cerr << "HepLorentzVector::plus() - "
      << "A zero vector used as reference to LorentzVector plus-part"
      << std::endl;
    return ee;
  }
  return ee + pp.dot(ref)/r;
} /* plus */

double HepLorentzVector::minus (const Hep3Vector & ref) const {
  double r = ref.mag();
  if (r == 0) {
    std::cerr << "HepLorentzVector::minus() - "
      << "A zero vector used as reference to LorentzVector minus-part"
      << std::endl;
    return ee;
  }
  return ee - pp.dot(ref)/r;
} /* plus */

HepLorentzVector HepLorentzVector::rest4Vector() const {
  return HepLorentzVector (0, 0, 0, (t() < 0.0 ? -m() : m()));
}

//-********
// beta
// gamma
//-********

double HepLorentzVector::beta() const {
  if (ee == 0) {
    if (pp.mag2() == 0) {
      return 0;
    } else {
      std::cerr << "HepLorentzVector::beta() - "
        << "beta computed for HepLorentzVector with t=0 -- infinite result"
        << std::endl;
      return 1./ee;
    }
  }
//  if (restMass2() <= 0) {
//    std::cerr << "HepLorentzVector::beta() - "
//      << "beta computed for a non-timelike HepLorentzVector" << std::endl;
//        // result will make analytic sense but is physically meaningless
//  }
  return std::sqrt (pp.mag2() / (ee*ee)) ;
} /* beta */

double HepLorentzVector::gamma() const {
  double v2 = pp.mag2();
  double t2 = ee*ee;
  if (ee == 0) {
    if (pp.mag2() == 0) {
      return 1;
    } else {
      std::cerr << "HepLorentzVector::gamma() - "
        << "gamma computed for HepLorentzVector with t=0 -- zero result"
        << std::endl;
      return 0;
    }
  }
  if (t2 < v2) {
    std::cerr << "HepLorentzVector::gamma() - "
      << "gamma computed for a spacelike HepLorentzVector -- imaginary result"
      << std::endl;
        // analytic result would be imaginary.
    return 0;
//  } else if ( t2 == v2 ) {
//    std::cerr << "HepLorentzVector::gamma() - "
//      << "gamma computed for a lightlike HepLorentzVector -- infinite result"
//      << std::endl;
  }
  return 1./std::sqrt(1. - v2/t2 );
} /* gamma */


//-***************
// rapidity
// pseudorapidity
// eta
//-***************

double HepLorentzVector::rapidity() const {
  double z1 = pp.getZ();
//  if (std::fabs(ee) == std::fabs(z1)) {
//    std::cerr << "HepLorentzVector::rapidity() - "
//      << "rapidity for 4-vector with |E| = |Pz| -- infinite result"
//      << std::endl;
//  }
  if (std::fabs(ee) < std::fabs(z1)) {
    std::cerr << "HepLorentzVector::rapidity() - "
      << "rapidity for spacelike 4-vector with |E| < |Pz| -- undefined"
      << std::endl;
    return 0;
  }
  double q = (ee + z1) / (ee - z1);
        //-| This cannot be negative now, since both numerator
        //-| and denominator have the same sign as ee.
  return .5 * std::log(q);
} /* rapidity */

double HepLorentzVector::rapidity(const Hep3Vector & ref) const {
  double r = ref.mag2();
  if (r == 0) {
    std::cerr << "HepLorentzVector::rapidity() - "
      << "A zero vector used as reference to LorentzVector rapidity"
      << std::endl;
    return 0;
  }
  double vdotu = pp.dot(ref)/std::sqrt(r);
//  if (std::fabs(ee) == std::fabs(vdotu)) {
//    std::cerr << "HepLorentzVector::rapidity() - "
//      << "rapidity for 4-vector with |E| = |Pu| -- infinite result"
//      << std::endl;
//  }
  if (std::fabs(ee) < std::fabs(vdotu)) {
    std::cerr << "HepLorentzVector::rapidity() - "
      << "rapidity for spacelike 4-vector with |E| < |P*ref| -- undefined "
      << std::endl;
    return 0;
  }
  double q = (ee + vdotu) / (ee - vdotu);
  return .5 * std::log(q);
} /* rapidity(ref) */

double HepLorentzVector::coLinearRapidity() const {
  double v1 = pp.mag();
//  if (std::fabs(ee) == std::fabs(v1)) {
//    std::cerr << "HepLorentzVector::coLinearRapidity() - "
//      << "co-Linear rapidity for 4-vector with |E| = |P| -- infinite result"
//      << std::endl;
//  }
  if (std::fabs(ee) < std::fabs(v1)) {
    std::cerr << "HepLorentzVector::coLinearRapidity() - "
      << "co-linear rapidity for spacelike 4-vector -- undefined"
      << std::endl;
    return 0;
  }
  double q = (ee + v1) / (ee - v1);
  return .5 * std::log(q);
} /* rapidity */

//-*************
// invariantMass
//-*************

double HepLorentzVector::invariantMass(const HepLorentzVector & w) const {
  double m1 = invariantMass2(w);
  if (m1 < 0) {
    // We should find out why:
    if ( ee * w.ee < 0 ) {
      std::cerr << "HepLorentzVector::invariantMass() - "
        << "invariant mass meaningless: \n"
        << "a negative-mass input led to spacelike 4-vector sum" << std::endl;
      return 0;
    } else if ( (isSpacelike() && !isLightlike()) ||
                (w.isSpacelike() && !w.isLightlike()) ) {
        std::cerr << "HepLorentzVector::invariantMass() - "
          << "invariant mass meaningless because of spacelike input"
          << std::endl;
      return 0;
    } else {
      // Invariant mass squared for a pair of timelike or lightlike vectors
      // mathematically cannot be negative.  If the vectors are within the
      // tolerance of being lightlike or timelike, we can assume that prior
      // or current roundoffs have caused the negative result, and return 0
      // without comment.
      return 0;
    }
  }
  return (ee+w.ee >=0 ) ? std::sqrt(m1) : - std::sqrt(m1);
} /* invariantMass */

//-***************
// findBoostToCM
//-***************

Hep3Vector HepLorentzVector::findBoostToCM() const {
  return -boostVector();
} /* boostToCM() */

Hep3Vector HepLorentzVector::findBoostToCM (const HepLorentzVector & w) const {
  double t1 = ee + w.ee;
  Hep3Vector v1 = pp + w.pp;
  if (t1 == 0) {
    if (v1.mag2() == 0) {
      return Hep3Vector(0,0,0);
    } else {
      std::cerr << "HepLorentzVector::findBoostToCM() - "
        << "boostToCM computed for two 4-vectors with combined t=0 -- "
        << "infinite result" << std::endl;
      return Hep3Vector(v1*(1./t1)); // Yup, 1/0 -- that is how we return infinity
    }
  }
//  if (t1*t1 - v1.mag2() <= 0) {
//    std::cerr << "HepLorentzVector::findBoostToCM() - "
//      << "boostToCM  computed for pair of HepLorentzVectors with non-timelike sum"
//      << std::endl;
//        // result will make analytic sense but is physically meaningless
//  }
  return Hep3Vector(v1 * (-1./t1));
} /* boostToCM(w) */

}  // namespace CLHEP

