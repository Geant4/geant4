#include <math.h>
#include "pair.h"
#include "vector"

namespace InuclSpecialFunctions {

double bindingEnergyExact(double A, double Z);

double bindingEnergyKummel(double A, double Z);

double bindingEnergy(double A, double Z);

double bindingEnergyAsymptotic(double A, double Z);

double FermiEnergy(double A, double Z, int ntype);
  
pair<vector<double>,vector<double> > paraMaker(double Z);

pair<double,double> paraMakerTruncated(double Z); 

double getAL(double A);
 
double csNN(double e);

double csPN(double e);

double inuclRndm();

double randomGauss(double sigma);

double randomPHI();

pair<double,double> randomCOS_SIN();

double nucleiLevelDensity(double a);

vector<double> generateWithFixedTheta(double ct, double p);

}
