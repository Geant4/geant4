#ifndef _HEATBATH_
#define _HEATBATH_

#ifdef IS_GCC
#pragma interface
#endif

#include "newvector.hh"
#include "Memory.hh"


template<class dist>
class HeatBath
{
public:
  HeatBath() : Gas(0) {}
  Vektor3 operator()(double T,double M);
private:
  dist* Gas;
};

template<class dist>
Vektor3 HeatBath<dist>::operator()(double T,double M)
{
  Gas = NEW dist(T,M);
  double r = rand_gen::Random();
  double plen = Gas->getValue();
  Vektor3 p = Vektor3::isotropy(plen);
  delete Gas;
  return p;
}

#endif
