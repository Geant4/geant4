#ifndef __MATHTOOLS__
#define __MATHTOOLS__

#include "double.hh"
#include "globals.hh"

class IntervalRange
{
public:
  IntervalRange() : a(0),b(0) {}
  IntervalRange(const IntervalRange& xx) : a(xx.a),b(xx.b) {}
  IntervalRange(REAL x0,REAL x1) : a(x0),b(x1) {}
  void set(REAL x0,REAL x1) { a = x0; b = x1; }
  IntervalRange& operator=(const IntervalRange& xx) 
       { set(xx.a,xx.b); return *this; }

  REAL a,b;
};

class GaussIntegration
{
public:
  REAL Integral(REAL);
protected:
  GaussIntegration(REAL a) : x0(a) {}
private:
  virtual REAL ToBeIntegrated(REAL) = 0;
  REAL x0;
};

class GaussLaguerre
{
protected:
  REAL Integral();
private:
  virtual REAL ToBeIntegrated(REAL) = 0;
};

class InverseFunction
{
public:
  static bool reducePrecision;
  REAL Inverse(REAL) const;
  void setRange(REAL a,REAL b) { x.set(a,b); }
  virtual int DivideInterval(IntervalRange,REAL,REAL&,int,double eps = Double::Epsilon) const;
protected:
  InverseFunction(REAL a,REAL b,double eps = Double::Epsilon) 
    : x(a,b),Epsilon(eps) {}
protected:
  virtual REAL ToBeInverted(REAL) const = 0;
  virtual REAL FindRoot(REAL,REAL,REAL,double) const;
  IntervalRange x;
  double Epsilon;
};

class InverseFunctionWithNewton : public InverseFunction
{
protected:
  InverseFunctionWithNewton(REAL x0,double eps=Double::Epsilon) 
    : InverseFunction(x0,x0,eps) {}
  int DivideInterval(IntervalRange,REAL,REAL&,int,double = Double::Epsilon) const { return 0; }
protected:
  virtual REAL Derivative(REAL x) const { return (ToBeInverted(x+Epsilon)-ToBeInverted(x))/Epsilon; }
  REAL FindRoot(REAL,REAL,REAL,double) const;
};

#endif
