//REQUIRES DistributionFunction.C
#ifndef __DistributionFunction__
#define __DistributionFunction__

#include "Random.hh"
#include "MathTools.hh"
#include "Fallible.hh"

class DistributionBase
{
public:
  virtual REAL getValue() const = 0;
  virtual REAL f(REAL x) const = 0;
};

class GaussDeviates : public DistributionBase
{
  Fallible<REAL>* store;
  REAL sqrt_2pi;
public:
  GaussDeviates() : store(new Fallible<REAL>),sqrt_2pi(sqrt(2.0*mathConstants::Pi)) {}
  REAL getValue() const;
  REAL f(REAL x) const { return exp(-0.5*x*x)/sqrt_2pi; }
};

class DistributionFunction : public DistributionBase
{
  int NMAX;
  REAL RangeMin,RangeMax;
public:
  REAL getValue() const;
  int getMaxRejections() const { return NMAX; }

  // f must be normalized to 1!!! (Int[f,0,1] = 1)
protected:
  DistributionFunction() : NMAX(0),RangeMin(0.0),RangeMax(1e+30) {}
  DistributionFunction(REAL a,REAL b) : NMAX(0),RangeMin(a),RangeMax(b) {}
  virtual REAL integratedMajorante(REAL x) const = 0;
private:
  virtual REAL Majorante(REAL x) const = 0;
  virtual REAL inverseIntegratedMajorante(REAL x) const = 0;
};

class SimpleRejection : public DistributionFunction
{
protected:
  SimpleRejection(REAL a,REAL b) : DistributionFunction(a,b) {}
private:
  REAL Majorante(REAL) const  { return 1.0; }
  REAL integratedMajorante(REAL x) const { return x; }
  REAL inverseIntegratedMajorante(REAL x) const { return x; }
 
};

class InvertedRejection : public DistributionFunction,
			  private InverseFunctionWithNewton
{
protected:
  InvertedRejection(REAL x0) 
    : DistributionFunction(),
      InverseFunctionWithNewton(x0) {}
  InvertedRejection(REAL a,REAL b,REAL x0) 
    : DistributionFunction(a,b),
      InverseFunctionWithNewton(x0) {}
protected:
  REAL Derivative(REAL x) const { return f(x); }
  REAL ToBeInverted(REAL x) const { return integratedMajorante(x); }
  REAL inverseIntegratedMajorante(REAL x) const { return InverseFunctionWithNewton::Inverse(x); }
};

#endif
