#ifndef __POTBASE__
#define __POTBASE__

#include "MathTools.hh"

class PotentialBase
{
public:
  virtual double V(double r,double s) const = 0;
  virtual double V_prime(double r,double s) const = 0;
  virtual double V_inv(double r,double s) const = 0;
  void print(G4std::ostream& o,double min,double max,int n);
};

class Potential : public PotentialBase,private InverseFunction
{
  double eps,s0;
  REAL ToBeInverted(REAL x) const { return V(x,s0); }
public:
  Potential(double a,double b,double dx = 1e-6) : InverseFunction(a,b),eps(dx) {}
  virtual double V_prime(double r,double s) const { return (V(r+eps,s)-V(r,s))/eps; }
  virtual double V_inv(double r,double s) const { (double&)s0 = s; return Inverse(r); }
};

#endif
