#include "ThermDist.hh"
#include "bessel.hh"
#include "Error.hh"

Relativistic::Relativistic(REAL t_,REAL m_) 
  : InvertedRejection(t_),T(t_),m(m_),norm(m*m*T*bessk(2,m/T)) {}

REAL Relativistic::f(REAL p) const 
{ 
  return p*p*exp(-sqrt(p*p+m*m)/T)/norm; 
}

REAL Relativistic::Majorante(REAL p) const 
{ 
  return p*p*exp(-p/T)/norm; 
}

REAL Relativistic::integratedMajorante(REAL x) const 
{ 
  REAL pT = x/T; 
  return 2.0*pow(T,3)*(1.0-exp(-pT)*(0.5*sqr(pT)+pT+1.0))/norm; 
}


REAL Relativistic::getValue() const
{
  REAL y;
  try {
    y = InvertedRejection::getValue();
  }
  catch (...) {
    y = -T*log(0.1*rand_gen::Random()*norm/(2.0*pow(T,3))); 
  }
  return y;
}



