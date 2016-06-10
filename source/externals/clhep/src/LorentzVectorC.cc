// -*- C++ -*-
// $Id:$
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the HepLorentzVector class:
// Those methods originating with ZOOM dealing with comparison (other than
// isSpaceLike, isLightlike, isTimelike, which are in the main part.)
//
// 11/29/05 mf in deltaR, replaced the direct subtraction 
// pp.phi() - w.getV().phi() with pp.deltaRPhi(w.getV()) which behaves 
// correctly across the 2pi boundary.

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/LorentzVector.h"

#include <cmath>

namespace CLHEP  {

//-***********
// Comparisons
//-***********

int HepLorentzVector::compare (const HepLorentzVector & w) const {
  if       ( ee > w.ee ) {
    return 1;
  } else if ( ee < w.ee ) {
    return -1;
  } else {
    return ( pp.compare(w.pp) );
  }
} /* Compare */

bool HepLorentzVector::operator > (const HepLorentzVector & w) const {
        return (compare(w)  > 0);
}
bool HepLorentzVector::operator < (const HepLorentzVector & w) const {
        return (compare(w)  < 0);
}
bool HepLorentzVector::operator>= (const HepLorentzVector & w) const {
        return (compare(w) >= 0);
}
bool HepLorentzVector::operator<= (const HepLorentzVector & w) const {
        return (compare(w) <= 0);
}

//-********
// isNear
// howNear
//-********

bool HepLorentzVector::isNear(const HepLorentzVector & w, 
						double epsilon) const {
  double limit = std::fabs(pp.dot(w.pp));
  limit += .25*((ee+w.ee)*(ee+w.ee));
  limit *= epsilon*epsilon;
  double delta = (pp - w.pp).mag2();
  delta +=  (ee-w.ee)*(ee-w.ee);
  return (delta <= limit );
} /* isNear() */

double HepLorentzVector::howNear(const HepLorentzVector & w) const {
  double wdw = std::fabs(pp.dot(w.pp)) + .25*((ee+w.ee)*(ee+w.ee));
  double delta = (pp - w.pp).mag2() + (ee-w.ee)*(ee-w.ee);
  if ( (wdw > 0) && (delta < wdw)  ) {
    return std::sqrt (delta/wdw);
  } else if ( (wdw == 0) && (delta == 0) ) {
    return 0;
  } else {
    return 1;
  }
} /* howNear() */

//-*********
// isNearCM
// howNearCM
//-*********

bool HepLorentzVector::isNearCM
			(const HepLorentzVector & w, double epsilon) const {

  double tTotal = (ee + w.ee);
  Hep3Vector vTotal (pp + w.pp);
  double vTotal2 = vTotal.mag2();

  if ( vTotal2 >= tTotal*tTotal ) {
    // Either one or both vectors are spacelike, or the dominant T components
    // are in opposite directions.  So boosting and testing makes no sense;
    // but we do consider two exactly equal vectors to be equal in any frame,
    // even if they are spacelike and can't be boosted to a CM frame.
    return (*this == w);
  }

  if ( vTotal2 == 0 ) {  // no boost needed!
    return (isNear(w, epsilon));
  }

  // Find the boost to the CM frame.  We know that the total vector is timelike.

  double tRecip = 1./tTotal;
  Hep3Vector bboost ( vTotal * (-tRecip) );

        //-| Note that you could do pp/t and not be terribly inefficient since
        //-| SpaceVector/t itself takes 1/t and multiplies.  The code here saves
        //-| a redundant check for t=0.

  // Boost both vectors.  Since we have the same boost, there is no need
  // to repeat the beta and gamma calculation; and there is no question
  // about beta >= 1.  That is why we don't just call w.boosted().

  double b2 = vTotal2*tRecip*tRecip;

  double ggamma = std::sqrt(1./(1.-b2));
  double boostDotV1 = bboost.dot(pp);
  double gm1_b2 = (ggamma-1)/b2;

  HepLorentzVector w1 ( pp   + ((gm1_b2)*boostDotV1+ggamma*ee) * bboost,
                     ggamma * (ee + boostDotV1) );

  double boostDotV2 = bboost.dot(w.pp);
  HepLorentzVector w2 ( w.pp + ((gm1_b2)*boostDotV2+ggamma*w.ee) * bboost,
                     ggamma * (w.ee + boostDotV2) );

  return (w1.isNear(w2, epsilon));

} /* isNearCM() */

double HepLorentzVector::howNearCM(const HepLorentzVector & w) const {

  double tTotal = (ee + w.ee);
  Hep3Vector vTotal (pp + w.pp);
  double vTotal2 = vTotal.mag2();

  if ( vTotal2 >= tTotal*tTotal ) {
    // Either one or both vectors are spacelike, or the dominant T components
    // are in opposite directions.  So boosting and testing makes no sense;
    // but we do consider two exactly equal vectors to be equal in any frame,
    // even if they are spacelike and can't be boosted to a CM frame.
    if (*this == w) {
      return 0;
    } else {
      return 1;
    }
  }

  if ( vTotal2 == 0 ) {  // no boost needed!
    return (howNear(w));
  }

  // Find the boost to the CM frame.  We know that the total vector is timelike.

  double tRecip = 1./tTotal;
  Hep3Vector bboost ( vTotal * (-tRecip) );

        //-| Note that you could do pp/t and not be terribly inefficient since
        //-| SpaceVector/t itself takes 1/t and multiplies.  The code here saves
        //-| a redundant check for t=0.

  // Boost both vectors.  Since we have the same boost, there is no need
  // to repeat the beta and gamma calculation; and there is no question
  // about beta >= 1.  That is why we don't just call w.boosted().

  double b2 = vTotal2*tRecip*tRecip;
//  if ( b2 >= 1 ) {			// NaN-proofing
//    std::cerr << "HepLorentzVector::howNearCM() - "
//	<< "boost vector in howNearCM appears to be tachyonic" << std::endl;
//  }
  double ggamma = std::sqrt(1./(1.-b2));
  double boostDotV1 = bboost.dot(pp);
  double gm1_b2 = (ggamma-1)/b2;

  HepLorentzVector w1 ( pp   + ((gm1_b2)*boostDotV1+ggamma*ee) * bboost,
                     ggamma * (ee + boostDotV1) );

  double boostDotV2 = bboost.dot(w.pp);
  HepLorentzVector w2 ( w.pp + ((gm1_b2)*boostDotV2+ggamma*w.ee) * bboost,
                     ggamma * (w.ee + boostDotV2) );

  return (w1.howNear(w2));

} /* howNearCM() */

//-************
// deltaR
// isParallel
// howParallel
// howLightlike
//-************

double HepLorentzVector::deltaR ( const HepLorentzVector & w ) const {

  double a = eta() - w.eta();
  double b = pp.deltaPhi(w.getV());

  return std::sqrt ( a*a + b*b );

} /* deltaR */

// If the difference (in the Euclidean norm) of the normalized (in Euclidean
// norm) 4-vectors is small, then those 4-vectors are considered nearly
// parallel.

bool HepLorentzVector::isParallel (const HepLorentzVector & w, double epsilon) const {
  double norm = euclideanNorm();
  double wnorm = w.euclideanNorm();
  if ( norm == 0 ) {
    if ( wnorm == 0 ) {
      return true;
    } else {
      return false;
    }
  }
  if ( wnorm == 0 ) {
    return false;
  }
  HepLorentzVector w1 = *this / norm;
  HepLorentzVector w2 = w / wnorm;
  return ( (w1-w2).euclideanNorm2() <= epsilon*epsilon );
} /* isParallel */


double HepLorentzVector::howParallel (const HepLorentzVector & w) const {

  double norm = euclideanNorm();
  double wnorm = w.euclideanNorm();
  if ( norm == 0 ) {
    if ( wnorm == 0 ) {
      return 0;
    } else {
      return 1;
    }
  }
  if ( wnorm == 0 ) {
    return 1;
  }

  HepLorentzVector w1 = *this / norm;
  HepLorentzVector w2 = w / wnorm;
  double x1 = (w1-w2).euclideanNorm();
  return (x1 < 1) ? x1 : 1;

} /* howParallel */

double HepLorentzVector::howLightlike() const {
  double m1 = std::fabs(restMass2());
  double twoT2 = 2*ee*ee;
  if (m1 < twoT2) {
    return m1/twoT2;
  } else {
    return 1;
  }
} /* HowLightlike */

}  // namespace CLHEP
