#include "PotentialBase.hh"

class Linear : public PotentialBase
{
  double kappa;
public:
  Linear(double k_) : kappa(k_) { cerr << "Potential: Linear\n"; }
  virtual double V(double r,double s) const { return s*kappa*r; }
  virtual double V_prime(double r,double s) const { return s*kappa; }
  virtual double V_inv(double E,double s) const { return E/kappa; }
};

class Cornell : public Potential
{
  double e,C,kappa;
public:
  Cornell() : Potential(0.01,10),e(0.25),C(0),kappa(1) { cerr << "Potential: Cornell\n"; }
  virtual double V(double r,double s) const { return s*(C-e/r+kappa*r); }
  virtual double V_prime(double r,double s) const { return s*(e/(r*r)+kappa); }
  virtual double V_inv(double E,double s) const { double x=(s*E-C)/(2.0*kappa); return x+sqrt(x*x+e/kappa); }
};

class InterQuark1 : public Potential
{
  REAL kappa,alpha,mu;
  REAL V1(REAL r) const { return -4.0/3.0*(alpha/r)*(1+alpha*3/(2*M_PI)*(11.0/3.0*log(mu*r)-11/3*0.5772+31.0/9.0))+kappa*r; }
  REAL V2(REAL r) const { return 4.0/3.0*(alpha/r)-kappa*r; }
public:
  InterQuark1() : Potential(0.01,10),kappa(0.9),mu(sqrt(kappa/0.197)),alpha(0.2) {}
  virtual double V(double r,double s) const { return (s>0) ? fabs(s)*V1(r) : fabs(s)*V2(r); }
  //  virtual double V_prime(double r,double s) const { return fabs(s)*4.0/3.0*alpha/(r*r)*(sign(s)-11*alpha/(2*M_PI)*(1-log(mu*r)))+s*kappa; }
};

class InterQuark : public Potential
{
  REAL kappa,alpha,mu,A,B,C,r0,V0;
  REAL V1(REAL r) const { return (r>=r0) ? -A/r*(B*log(mu*r)+C) : V0; }
  REAL V2(REAL r) const { return A/r; }
public:
  InterQuark(double k = 0.9) : Potential(0.1,10),kappa(k),mu(sqrt(kappa/0.197)),
    alpha(0.2),A(mathConstants::hc*4.0/3.0*alpha),B(alpha*3/(2*M_PI)*(11.0/3.0)),
    C(1+alpha*3/(2*M_PI)*(-11/3*0.5772+31.0/9.0)),
    r0(exp(1.0-C/B)/mu),V0(0) { V0 = V1(r0); cerr << "Potential: interquark\n"; }
  virtual double V(double r,double s) const { return ((s>0) ? s*V1(fabs(r)) : -s*V2(fabs(r)))+s*kappa*r; }
};
