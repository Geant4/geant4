#include "InuclSpecialFunctions.h"

double InuclSpecialFunctions::getAL(double A) {
  return 0.76 + 2.2/pow(A,0.333333);
}

double InuclSpecialFunctions::csNN(double e) {
  double snn;
  if(e < 40.) {
    snn = -1174.8/(e*e) + 3088.5/e + 5.3107;
  }
   else {
    snn = 93074./(e*e) - 11.148/e + 22.429;
  };
  return snn; 
}

double InuclSpecialFunctions::csPN(double e) {
  double spn;
  if(e < 40.) {
    spn = -5057.4/(e*e) + 9069.2/e + 6.9466;
  }
   else {
    spn = 239380./(e*e) + 1802./e + 27.147;
  };
  return spn; 
}

double InuclSpecialFunctions::FermiEnergy(double A, double Z, int ntype) {
// calculates the nuclei Fermi energy for 0 - neutron and 1 - proton

  const double C = 55.4;
  double Ef;
  if(ntype == 0) {
    Ef = C*pow((A - Z)/A, 0.666667);
  }
   else {
    Ef = C*pow(Z/A, 0.666667);
  };
  return Ef; 
}

double InuclSpecialFunctions::inuclRndm() { 
// PUT PROPER RNDM HERE
return drand48(); 
} 

double InuclSpecialFunctions::randomGauss(double sigma) {

const double eps = 1.e-6;
const double twopi = 6.2831854;

double r1 = inuclRndm();
r1 = r1 > eps ? r1 : eps;
double r2 = inuclRndm();
r2 = r2 > eps ? r2 : eps;
r2 = r2 < 1. - eps ? r2 : 1. - eps; 
return sigma*sin(twopi*r1)*sqrt(-2.*log(r2)); 

} 

double InuclSpecialFunctions::randomPHI() { 
  const double twopi = 6.2831853;
  return twopi*inuclRndm();
} 

pair<double,double> InuclSpecialFunctions::randomCOS_SIN() {
  double CT = 1. - 2.*inuclRndm();
  return pair<double,double>(CT,sqrt(1. - CT*CT));
}

vector<double> InuclSpecialFunctions::generateWithFixedTheta(double ct,
               double p) {

vector<double> momr(4);
double phi = randomPHI();
double pt = p*sqrt(fabs(1. - ct*ct));
vector<double> mom1(4);
momr[1] = pt*cos(phi);
momr[2] = pt*sin(phi);
momr[3] = p*ct;

return momr;

}
