#include "InuclSpecialFunctions.h"

double InuclSpecialFunctions::bindingEnergyAsymptotic(double A, double Z) {

// calculates the nuclei binding energy using smooth liquid high energy formula
  
  double X = (1.-2.*Z/A)*(1.-2.*Z/A);
  double X1 = pow(A,0.3333333);
  double X2 = X1*X1;
  double X3 = 1./X1;
  double X4 = 1./X2;
  double DM = 17.035*(1.-1.846*X) * A -
         25.8357*(1.-1.712*X) * X2 * ((1.-0.62025*X4)*(1.-0.62025*X4)) -
         0.779*Z*(Z-1.)*X3*(1.-1.5849*X4+1.2273/A+1.5772*X4*X4) +
         0.4328*pow(Z,1.333333) * X3 * (1.-0.57811*X3-0.14518*X4+0.496/A);
  return DM;

}
