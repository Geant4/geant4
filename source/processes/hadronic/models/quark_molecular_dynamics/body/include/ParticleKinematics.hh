#ifndef __PARTICLEKINEMATICS__
#define __PARTICLEKINEMATICS__

#include "newvector.hh"

double CMmomentum(double M,double m1,double m2);

Vektor4 Lorentz(const Vektor3& beta,const Vektor4& P);

#endif
