#ifndef __KINEMATICS__
#define __KINEMATICS__

#include "ParticleBase.hh"
#include "DistributionFunction.hh"
#include "Fallible.hh"
#include "Error.hh"

class FragmentationFunction : public SimpleRejection
{
  static REAL a,b;
  Fallible<REAL> mt2;
public:
  FragmentationFunction(REAL min_,REAL max_) 
    : SimpleRejection(min_,max_),mt2() {}
  FragmentationFunction(REAL min_,REAL max_,REAL mt2_) 
    : SimpleRejection(min_,max_),mt2(mt2_) {}
  void set_mt(REAL mt2_) { mt2 = mt2_; }
  REAL f(REAL z) const { return pow(1-z,a)*exp(-b*mt2/z); }
};

class Banerjee : public DistributionBase
{
  REAL a,b,ex;
public:
  Banerjee(REAL a_,REAL b_) : ex(exp(-0.5*zeta)) 
    {
      a = Int_f(a_);
      b = Int_f(b_);
    }
  REAL getValue() const { REAL px=rand_gen::Random(a,b); return 1.0-sqrt(-2.0/zeta*log(px+(1-px)*ex)); }
  static REAL zeta;
protected:
  REAL f(REAL x) const { return zeta/(1.0-ex)*(1.0-x)*exp(-0.5*zeta*sqr(1-x)); }
  REAL Int_f(REAL x) { return (exp(-0.5*zeta*sqr(1-x))-ex)/(1-ex); }
};

class Kinematics
{
  REAL pq,pq2,sigma,sigma_l;
  REAL Mt,Mt2,mQt,mQt2,eps,eps_min,eps_max_0,eps_max_1;
  const Vektor4& p0,&x0;
  Vektor4& x1,&x2,&p1,&p2;
  double m0,m1,m2;
public:
  typedef void (Kinematics::*chooseKin)(REAL&,REAL&);

  Kinematics(const Vektor4& x0_,Vektor4& x1_,Vektor4& x2_,
	     const Vektor4& p0_,Vektor4& p1_,Vektor4& p2_,
	     double m0_,double m1_,double m2_)
    : p0(p0_),p1(p1_),p2(p2_),x0(x0_),x1(x1_),x2(x2_),m0(m0_),m1(m1_),m2(m2_) {}
  void calculate(double);

  void kin1(REAL& z,REAL& l);
  void kin2(REAL& z,REAL& l);
  void kin3(REAL& z,REAL& l);
  
  static chooseKin whichKin;

  class NotPossible : public Error
  {
    int flavour;
  public:
    NotPossible() : flavour(0) {}
    NotPossible(int f) : flavour(f) {}
    operator int() const { return flavour; }
    void writeMessage(ostream& o) const { o << "Collision not possible..."; }
  };

private:
  REAL get_eps(REAL);
  REAL inv_eps(REAL,signed int);
};

#endif
