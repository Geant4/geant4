#include "DistributionFunction.hh"
#include "G4ios.hh"
#include <algo.h>

REAL GaussDeviates::getValue() const
{
  REAL y = 0;
  if ( store->isValid() ) {
    y = *store;
    store->invalidate();
  }
  else {
    REAL x1;
    REAL x2;
    REAL xx;
    do {
      x1 = rand_gen::Random(-1.0,1.0);
      x2 = rand_gen::Random(-1.0,1.0);
      xx = sqr(x1)+sqr(x2);
    }
    while ( xx > 1.0 || xx == 0.0 );
    REAL fac = sqrt(-2.0*log(xx)/xx);
    y = x1*fac;
    store->validate(x2*fac);
  }
  return y;
}

REAL DistributionFunction::getValue() const 
{
  int n = 0,counter=1000;
  REAL rmin = integratedMajorante(RangeMin);
  REAL rmax = integratedMajorante(RangeMax);
  REAL z;
  do {
    do {
      REAL r = rand_gen::Random(rmin,rmax);
      z = inverseIntegratedMajorante(r);
      //      G4cerr << z << "  " << f(z) << "  " << Majorante(z) << G4endl;
    }
    while ( z<RangeMin || z>RangeMax );
    ++n;
    //    G4cerr << z << "  " << Majorante(z) << "  " << f(z) << G4endl; 
  }
  while (Majorante(z)*rand_gen::Random() > f(z) && --counter);
  if ( !counter ) {
    G4cerr << "DistributionFunction::CounterOverflow...\n";
    throw "Counter overflow";
  }
  return z;
}

