#ifndef __THERMDIST__
#define __THERMDIST__

#include <math.h>
#include "DistributionFunction.hh"
#include "bessel.hh"

class Relativistic : public InvertedRejection
{
  REAL T,m,norm;
public:
  Relativistic(REAL t_,REAL m_);
  REAL getValue() const;
private:
  REAL f(REAL p) const;
  REAL Majorante(REAL p) const;
  REAL integratedMajorante(REAL x) const;
  //  REAL inverseIntegratedMajorante(REAL x) const;
};

#endif
