#include "Average.hh"

void  AverageVector::addEntry(const Vektor& x) {
  avg += x; 
  for (int i=1; i<=dim(avg2); i++)
    avg2[i] += x[i]*x[i]; 
  ++N; 
}

Vektor AverageVector::getError() const 
{ 
  Vektor err(dim(avg));
  if ( N>1 ) {
    Vektor x = avg/N; 
    for (int i=1; i<=dim(avg); i++) 
      err[i] = sqrt((avg2[i]-N*x[i]*x[i])/(N-1)); 
  } 
  return err; 
}

